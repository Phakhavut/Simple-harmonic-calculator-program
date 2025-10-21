from flask import Flask, request, jsonify, send_from_directory
import math
import copy
import re

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
            return val  # คืน str ถ้าไม่ใช่ float
    return val

def safe_sqrt(val):
    if val is None or val < 0:
        return None
    return math.sqrt(val)

def safe_div(num, den):
    if num is None or den is None or den == 0:
        return None
    return num / den

def relationship_engine(s1, s2, system="spring"):
    s2 = copy.deepcopy(s2)
    if system == "pendulum":
        if s1.get("L") and s1.get("T"):
            ratio = (s2.get("T", s1["T"]) / s1["T"]) ** 2
            g1 = s1.get("g", 9.81)
            g2 = s2.get("g", g1)
            ratio *= g1 / g2
            if s2.get("T") and not s2.get("L"):
                s2["L"] = s1["L"] * ratio
            elif s2.get("L") and not s2.get("T"):
                s2["T"] = s1["T"] * math.sqrt(s2["L"] / s1["L"] * (g2 / g1))
            elif s2.get("g") and not s2.get("T"):
                s2["T"] = 2 * math.pi * math.sqrt(s2["L"] / s2["g"])

    elif system == "spring":
        const = None
        if s1.get("m") and s1.get("k"):
            if s1.get("T"):
                const = s1["T"]**2 * s1["k"] / s1["m"]
            elif s1.get("omega"):
                const = (2 * math.pi / s1["omega"])**2 * s1["k"] / s1["m"]
            elif s1.get("f"):
                const = (1 / s1["f"])**2 * s1["k"] / s1["m"]
        if const is not None:
            if s2.get("m") and s2.get("T") and s2.get("k") is None:
                s2["k"] = const * s2["m"] / s2["T"]**2
            elif s2.get("k") and s2.get("T") and s2.get("m") is None:
                s2["m"] = s2["T"]**2 * s2["k"] / const
            elif s2.get("m") and s2.get("k") and s2.get("T") is None:
                s2["T"] = math.sqrt(const * s2["m"] / s2["k"])
    return s2

def parse_expression(expr, prev_context):
    if isinstance(expr, (int, float)): return expr
    if not isinstance(expr, str): return None
    expr = expr.strip()
    if expr == "=": return "="
    # แทนค่าตัวแปรจากสภาวะก่อนหน้า
    for k,v in prev_context.items():
        if v is not None:
            expr = re.sub(rf"\b{k}\b", str(v), expr)
    try:
        return eval(expr, {"__builtins__": {}}, {"math": math})
    except:
        return None

def finalize_state(s, system_type):
    max_iter = 20
    for _ in range(max_iter):
        old_s = copy.deepcopy(s)
        # คำนวณ omega
        if s.get("omega") is None:
            if s.get("T") is not None and s.get("T") != 0:
                s["omega"] = safe_div(2 * math.pi, s["T"])
            elif s.get("f") is not None:
                s["omega"] = 2 * math.pi * s["f"]
            elif system_type == 1 and s.get("k") is not None and s.get("m") is not None and s.get("m") != 0:
                s["omega"] = safe_sqrt(safe_div(s["k"], s["m"]))
            elif system_type == 2 and s.get("g") is not None and s.get("L") is not None and s.get("L") != 0:
                s["omega"] = safe_sqrt(safe_div(s["g"], s["L"]))
        # T
        if s.get("T") is None:
            if s.get("omega") is not None and s.get("omega") != 0:
                s["T"] = safe_div(2 * math.pi, s["omega"])
            elif s.get("f") is not None and s.get("f") != 0:
                s["T"] = safe_div(1, s["f"])
        # f
        if s.get("f") is None:
            if s.get("omega") is not None and s.get("omega") != 0:
                s["f"] = safe_div(s["omega"], 2 * math.pi)
            elif s.get("T") is not None and s.get("T") != 0:
                s["f"] = safe_div(1, s["T"])
        # sigmaF
        if s.get("sigmaF") is None:
            if system_type == 1 and s.get("k") is not None and s.get("x") is not None:
                s["sigmaF"] = s["k"] * s["x"]
            elif s.get("m") is not None and s.get("g") is not None:
                s["sigmaF"] = s["m"] * s["g"]
        # คำนวณ k/m สำหรับ spring
        if system_type == 1:
            if s.get("k") is None and s.get("m") is not None and s.get("omega") is not None:
                s["k"] = s["m"] * s["omega"]**2
            if s.get("m") is None and s.get("k") is not None and s.get("omega") is not None and s.get("omega") != 0:
                s["m"] = safe_div(s["k"], s["omega"]**2)
        # A
        if s.get("A") is None:
            if s.get("v") is not None and s.get("x") is not None and s.get("omega") is not None and s.get("omega") != 0:
                val = safe_div(s["v"]**2, s["omega"]**2)
                if val is not None:
                    s["A"] = safe_sqrt(val + s["x"]**2)
            elif s.get("Vmax") is not None and s.get("omega") is not None and s.get("omega") != 0:
                s["A"] = safe_div(s["Vmax"], s["omega"])
            elif s.get("a_max") is not None and s.get("omega") is not None and s.get("omega") != 0:
                s["A"] = safe_div(s["a_max"], s["omega"]**2)
            elif s.get("E") is not None and s.get("k") is not None and s.get("k") != 0:
                s["A"] = safe_sqrt(safe_div(2 * s["E"], s["k"]))
        # v
        if s.get("v") is None:
            if s.get("omega") is not None and s.get("A") is not None and s.get("x") is not None:
                val = s["A"]**2 - s["x"]**2
                if val >= 0:
                    s["v"] = s["omega"] * safe_sqrt(val)
        # Vmax
        if s.get("Vmax") is None:
            if s.get("omega") is not None and s.get("A") is not None:
                s["Vmax"] = s["omega"] * s["A"]
            elif s.get("E") is not None and s.get("m") is not None and s.get("m") != 0:
                s["Vmax"] = safe_sqrt(safe_div(2 * s["E"], s["m"]))
        # a_max
        if s.get("a_max") is None:
            if s.get("omega") is not None and s.get("A") is not None:
                s["a_max"] = s["omega"]**2 * s["A"]
            elif s.get("Vmax") is not None and s.get("omega") is not None:
                s["a_max"] = s["Vmax"] * s["omega"]
        # a
        if s.get("a") is None and s.get("omega") is not None and s.get("x") is not None:
            s["a"] = - s["omega"]**2 * s["x"]
        # PE
        if s.get("PE") is None and s.get("k") is not None and s.get("x") is not None:
            s["PE"] = 0.5 * s["k"] * s["x"]**2
        # KE
        if s.get("KE") is None:
            if s.get("m") is not None and s.get("v") is not None:
                s["KE"] = 0.5 * s["m"] * s["v"]**2
            elif s.get("E") is not None and s.get("PE") is not None:
                s["KE"] = s["E"] - s["PE"]
        # E
        if s.get("E") is None:
            if s.get("k") is not None and s.get("A") is not None:
                s["E"] = 0.5 * s["k"] * s["A"]**2
            elif s.get("m") is not None and s.get("Vmax") is not None:
                s["E"] = 0.5 * s["m"] * s["Vmax"]**2
            elif s.get("KE") is not None and s.get("PE") is not None:
                s["E"] = s["KE"] + s["PE"]
        # PE from E - KE
        if s.get("PE") is None and s.get("E") is not None and s.get("KE") is not None:
            s["PE"] = s["E"] - s["KE"]
        if old_s == s:
            break
    return s

def process_states(raw_states, system="spring"):
    system_type = 1 if system.lower().startswith("s") else 2
    states = []
    prev_values = {}
    s_prev_fin = None

    for i, raw in enumerate(raw_states):
        s = {}
        for k in NUMERIC_KEYS:
            val = to_num(raw.get(k))
            if isinstance(val, str) and val != "=":
                val = parse_expression(val, prev_values)
            if val == "=":
                val = prev_values.get(k)
            s[k] = val
        if i > 0 and s_prev_fin is not None:
            s = relationship_engine(s_prev_fin, s, system)
        if system_type == 2 and s.get("g") is None:
            s["g"] = 9.81
        s_fin = finalize_state(s, system_type)
        states.append(s_fin)
        s_prev_fin = s_fin
        for k, v in s_fin.items():
            if v is not None:
                prev_values[k] = v
    results = []
    for s in states:
        out = {k: s.get(k) for k in NUMERIC_KEYS if s.get(k) is not None}  # Filter None เพื่อ clean output
        out["warnings"] = []  # สำหรับ JS
        out["explanation"] = ""  # สำหรับ JS
        results.append(out)
    return results


@app.route("/")
def serve_html():
    return send_from_directory(".", "index.html")

@app.route("/solve", methods=["POST"])
def solve():
    data = request.get_json()
    system = data.get("system", "spring")
    raw_states = data.get("states", [])
    try:
        results = process_states(raw_states, system)
        return jsonify({"system": system, "results": results})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)
