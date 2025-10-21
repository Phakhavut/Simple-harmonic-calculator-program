# shm_solver_v3.py
import math
import copy
import re

NUMERIC_KEYS = [
    "m","k","L","g","sigmaF","x","v","Vmax","a_max","a",
    "KE","PE","E","A","omega","f","T"
]

def to_num(val):
    """แปลงค่า string เป็น float หรือ '=' หรือ None"""
    if val is None:
        return None
    if isinstance(val, (int, float)):
        return val
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
    return None


def safe_sqrt(val):
    if val is None or val < 0:
        return None
    return math.sqrt(val)


def safe_div(num, den):
    if num is None or den is None or den == 0:
        return None
    return num / den


def eval_param(expr, prev_values):
    """รองรับการใช้ตัวแปรเช่น 2x, 3*T, L/2"""
    if expr is None:
        return None
    if isinstance(expr, (int, float)):
        return expr
    expr = expr.strip()
    if expr == "":
        return None
    for key, val in prev_values.items():
        if val is not None:
            expr = re.sub(rf"\b{key}\b", str(val), expr)
    try:
        return float(eval(expr, {"__builtins__": {}}))
    except:
        return None


def finalize_state(s, system_type):
    """คำนวณค่าต่างๆของแต่ละสภาวะ"""
    for _ in range(8):
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

        # k/m (spring)
        if system_type == 1:
            if s.get("k") is None and s.get("m") not in (None,0) and s.get("omega") not in (None,0):
                s["k"] = s["m"] * s["omega"]**2
            if s.get("m") is None and s.get("k") not in (None,0) and s.get("omega") not in (None,0):
                s["m"] = safe_div(s["k"], s["omega"]**2)

        # A
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
            elif s.get("E") is not None and s.get("m") not in (None,0) and s.get("omega") not in (None,0):
                s["A"] = safe_sqrt(2*s["E"]/(s["m"]*s["omega"]**2))

        # v
        if s.get("v") is None:
            if s.get("omega") not in (None,0) and s.get("A") is not None and s.get("x") is not None:
                val = s["omega"]**2 * (s["A"]**2 - s["x"]**2)
                if val >= 0:
                    s["v"] = safe_sqrt(val)
            elif s.get("KE") is not None and s.get("m") not in (None,0):
                s["v"] = safe_sqrt(2*s["KE"]/s["m"])
            elif s.get("E") is not None and s.get("PE") is not None and s.get("m") not in (None,0):
                ke = s["E"] - s["PE"]
                if ke >= 0:
                    s["v"] = safe_sqrt(2*ke/s["m"])

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
    """ประมวลผลหลายสภาวะ + การแปรผัน + ตัวแปรคงที่ร่วม"""
    system_type = 1 if system.lower().startswith("s") else 2
    states = []
    prev_values = {}

    for raw in raw_states:
        s = {}
        for k in NUMERIC_KEYS:
            raw_val = raw.get(k)
            val = to_num(raw_val)
            if val == "=":
                val = prev_values.get(k)
            elif val is None and isinstance(raw_val, str):
                val = eval_param(raw_val, prev_values)
            elif len(states) >= 1 and val is None:
                # ถ้าไม่ได้กำหนดและมี state ก่อนหน้า → ใช้ค่าเดิม
                val = prev_values.get(k)
            s[k] = val
            if val is not None:
                prev_values[k] = val
        states.append(s)

    results = []
    for i, s in enumerate(states):
        s_copy = copy.deepcopy(s)
        if system_type == 2 and s_copy.get("g") is None:
            s_copy["g"] = 9.81

        if i > 0:
            prev = results[-1]
            # 🔸 ลูกตุ้ม: T² ∝ L
            if system_type == 2:
                if s_copy.get("L") is None and s_copy.get("T") is not None:
                    s_copy["L"] = prev["L"] * (s_copy["T"]/prev["T"])**2
                elif s_copy.get("T") is None and s_copy.get("L") is not None:
                    s_copy["T"] = prev["T"] * math.sqrt(s_copy["L"]/prev["L"])
            # 🔸 สปริง: T² ∝ m/k
            elif system_type == 1:
                if s_copy.get("T") is None and s_copy.get("m") is not None and s_copy.get("k") is not None:
                    s_copy["T"] = 2*math.pi*safe_sqrt(s_copy["m"]/s_copy["k"])
                elif s_copy.get("k") is None and s_copy.get("m") is not None and s_copy.get("T") is not None:
                    s_copy["k"] = (4*math.pi**2*s_copy["m"])/(s_copy["T"]**2)
                elif s_copy.get("m") is None and s_copy.get("k") is not None and s_copy.get("T") is not None:
                    s_copy["m"] = s_copy["k"] * (s_copy["T"]**2)/(4*math.pi**2)

        s_fin = finalize_state(s_copy, system_type)
        out = {k: s_fin.get(k) for k in NUMERIC_KEYS}
        results.append(out)

    return results


# 🧪 ตัวอย่างทดสอบ
if __name__ == "__main__":
    import json

    print("---- ลูกตุ้ม ----")
    pendulum = [
        {"L": 0.25, "T": 1.0},
        {"T": 3.0}  # ไม่ใส่ L → โปรแกรมคำนวณ L=2.25
    ]
    print(json.dumps(process_states(pendulum, "pendulum"), indent=2))

    print("\n---- สปริง ----")
    spring = [
        {"m": 0.5, "k": 200, "T": 1.0},
        {"T": 2.0}  # ไม่ใส่ k,m → โปรแกรมใช้ค่าเดิมและตรวจการแปรผัน
    ]
    print(json.dumps(process_states(spring, "spring"), indent=2))
