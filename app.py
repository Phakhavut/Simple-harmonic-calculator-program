from flask import Flask, request, jsonify, send_from_directory
import math
import copy

app = Flask(__name__)

NUMERIC_KEYS = [
    "m","k","L","g","sigmaF","x","v","Vmax","a_max","a",
    "KE","PE","E","A","omega","f","T"
]

def to_num(val):
    if val is None:
        return None
    if isinstance(val, str):
        val = val.strip()
        if val == "":
            return None
        if val == "=":
            return "="
        try:
            return float(val)
        except:
            return None
    return val

def safe_sqrt(val):
    if val is None or val < 0:
        return None
    return math.sqrt(val)

def safe_div(num, den):
    if num is None or den is None or den == 0:
        return None
    return num / den

def finalize_state(s, system_type):
    for _ in range(8):
        # คำนวณ omega
        if s.get("omega") is None:
            if s.get("T") not in (None, 0):
                s["omega"] = 2 * math.pi / s["T"]
            elif s.get("f") not in (None, 0):
                s["omega"] = 2 * math.pi * s["f"]
            elif system_type == 1 and s.get("k") not in (None, 0) and s.get("m") not in (None, 0):
                s["omega"] = safe_sqrt(s["k"]/s["m"])
            elif system_type == 2 and s.get("g") not in (None, 0) and s.get("L") not in (None, 0):
                s["omega"] = safe_sqrt(s["g"]/s["L"])
        # sigmaF
        if s.get("sigmaF") is None:
            if system_type == 1 and s.get("k") is not None and s.get("x") is not None:
                s["sigmaF"] = s["k"] * s["x"]
            elif s.get("m") is not None and s.get("g") is not None:
                s["sigmaF"] = s["m"] * s["g"]
        # คำนวณค่าต่าง ๆ ตามสูตร
        if system_type == 1:
            if s.get("k") is None and s.get("m") not in (None,0) and s.get("omega") not in (None,0):
                s["k"] = s["m"] * s["omega"]**2
            if s.get("m") is None and s.get("k") not in (None,0) and s.get("omega") not in (None,0):
                s["m"] = safe_div(s["k"], s["omega"]**2)
        if s.get("A") is None:
            if s.get("v") is not None and s.get("x") is not None and s.get("omega") not in (None,0):
                val = safe_div(s["v"]**2, s["omega"]**2)
                if val is not None:
                    s["A"] = safe_sqrt(val + s["x"]**2)
            elif s.get("Vmax") is not None and s.get("omega") not in (None,0):
                s["A"] = safe_div(s["Vmax"], s["omega"])
            elif s.get("a_max") is not None and s.get("omega") not in (None,0):
                s["A"] = safe_div(s["a_max"], s["omega"]**2)
            elif s.get("E") is not None and s.get("k") not in (None,0):
                s["A"] = safe_sqrt(2*s["E"]/s["k"])
        if s.get("v") is None:
            if s.get("omega") not in (None,0) and s.get("A") is not None and s.get("x") is not None:
                val = s["omega"]**2 * (s["A"]**2 - s["x"]**2)
                if val >= 0:
                    s["v"] = safe_sqrt(val)
        if s.get("Vmax") is None and s.get("omega") not in (None,0) and s.get("A") is not None:
            s["Vmax"] = s["omega"]*s["A"]
        if s.get("a_max") is None and s.get("omega") not in (None,0) and s.get("A") is not None:
            s["a_max"] = s["omega"]**2*s["A"]
        if s.get("a") is None and s.get("omega") not in (None,0) and s.get("x") is not None:
            s["a"] = - s["omega"]**2 * s["x"]
        if s.get("PE") is None and s.get("k") not in (None,0) and s.get("x") is not None:
            s["PE"] = 0.5*s["k"]*s["x"]**2
        if s.get("E") is None and s.get("k") not in (None,0) and s.get("A") is not None:
            s["E"] = 0.5*s["k"]*s["A"]**2
        if s.get("KE") is None and s.get("m") not in (None,0) and s.get("v") is not None:
            s["KE"] = 0.5*s["m"]*s["v"]**2
        if s.get("KE") is None and s.get("E") is not None and s.get("PE") is not None:
            s["KE"] = s["E"] - s["PE"]
        if s.get("E") is None and s.get("KE") is not None and s.get("PE") is not None:
            s["E"] = s["KE"] + s["PE"]
        if s.get("f") is None and s.get("omega") not in (None,0):
            s["f"] = s["omega"]/(2*math.pi)
        if s.get("T") is None and s.get("f") not in (None,0):
            s["T"] = 1/s["f"]
    return s

def process_states(raw_states, system="spring"):
    system_type = 1 if system.lower().startswith("s") else 2
    states = []
    prev_values = {}

    for raw in raw_states:
        s = {}
        for k in NUMERIC_KEYS:
            val = to_num(raw.get(k))
            if val == "=":
                val = prev_values.get(k)
            s[k] = val
            if val is not None:
                prev_values[k] = val
        states.append(s)

    results = []
    for s in states:
        s_copy = copy.deepcopy(s)
        if system_type == 2 and s_copy.get("g") is None:
            s_copy["g"] = 9.81
        s_fin = finalize_state(s_copy, system_type)
        out = {k: s_fin.get(k) for k in NUMERIC_KEYS}
        results.append(out)
    return results


@app.route("/")
def serve_html():
    return send_from_directory(".", "index.html")

@app.route("/solve", methods=["POST"])  # ตรงกับ fetch('/solve')
def solve():
    data = request.get_json()
    system = data.get("system", "spring")
    raw_states = data.get("states", [])
    results = process_states(raw_states, system)
    return jsonify({"system": system, "results": results})

if __name__ == "__main__":
    app.run(debug=True)
