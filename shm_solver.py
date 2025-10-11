# shm_solver.py
import math
import copy

NUMERIC_KEYS = [
    "m","k","L","g","sigmaF","x","v","Vmax","a_max","a",
    "KE","PE","E","A","omega","f","T"
]

def is_equal_token(val):
    return isinstance(val, str) and val.strip() == "="

def to_num(val):
    # Accept numeric or None or "=" string
    if val is None or (isinstance(val, str) and val.strip() == ""):
        return None
    if is_equal_token(val):
        return "="
    try:
        return float(val)
    except:
        return None

def fmt(v):
    if v is None:
        return None
    try:
        return float(v)
    except:
        return v

def warn_if_inconsistent(s, system_type):
    warns = []
    # 1) ω vs √(k/m) for spring
    if system_type == 1 and s.get("k") and s.get("m") and s.get("omega"):
        if s["m"] != 0:
            expected = math.sqrt(s["k"]/s["m"])
            if abs(s["omega"] - expected) > 0.01 * max(expected, 1e-9):
                warns.append(f"ω ({s['omega']:.6g}) inconsistent with √(k/m) ({expected:.6g})")
    # 2) f vs ω/(2π)
    if s.get("omega") and s.get("f"):
        expected = s["omega"]/(2*math.pi)
        if abs(s["f"] - expected) > 0.01 * max(expected, 1e-9):
            warns.append(f"f ({s['f']:.6g}) inconsistent with ω/(2π) ({expected:.6g})")
    # 3) T vs 1/f
    if s.get("f") and s.get("T"):
        expected = 1.0 / s["f"] if s["f"] != 0 else None
        if expected and abs(s["T"] - expected) > 0.01 * max(expected, 1e-9):
            warns.append(f"T ({s['T']:.6g}) inconsistent with 1/f ({expected:.6g})")
    # 4) kx vs sigmaF if both present
    if s.get("k") and s.get("x") and s.get("sigmaF"):
        expected = s["k"] * s["x"]
        if abs(s["sigmaF"] - expected) > 0.01 * max(abs(expected),1e-9):
            warns.append(f"ΣF ({s['sigmaF']:.6g}) inconsistent with k·x ({expected:.6g})")
    # 5) energy/kinetic relations (rough)
    if s.get("m") and s.get("v") and s.get("KE"):
        expected = 0.5 * s["m"] * s["v"]**2
        if abs(s["KE"] - expected) > 0.01 * max(expected,1e-9):
            warns.append(f"KE ({s['KE']:.6g}) inconsistent with ½ m v² ({expected:.6g})")
    return warns

def explain_physics(system_type):
    if system_type == 1:
        txt = (
            "สปริง (mass-spring)\n"
            "- ω = √(k/m)\n"
            "- T = 2π√(m/k)\n"
            "- f = 1/T = ω/(2π)\n"
            "- ΣF (at equilibrium) = k x = m g\n"
            "- x(t) = A cos(ωt + φ) หรือ A sin(...)\n"
            "- v = ± ω √(A² - x²)\n"
            "- a = -ω² x\n"
            "- Vmax = ω A, a_max = ω² A\n"
            "- KE = ½ m v², PE = ½ k x², E = ½ k A²\n"
        )
    else:
        txt = (
            "ลูกตุ้ม (simple pendulum, small angle)\n"
            "- ω = √(g/L)\n"
            "- T = 2π √(L/g)\n"
            "- f = 1/T\n"
            "- x (linear small-angle) ~ A sin(ωt + φ) (arc-length approximation)\n"
            "- v = ± ω √(A² - x²), a = -ω² x (for SHM approx)\n"
        )
    return txt

# core inference & finalize
def infer_global_params(states, system_type):
    # If pendulum, infer g from any T/L or omega/L if present
    if system_type == 2:
        g_vals = []
        for s in states:
            if s.get("g") is not None:
                g_vals.append(s["g"])
            if s.get("T") and s.get("L"):
                try:
                    g_vals.append(4 * math.pi**2 * s["L"] / (s["T"]**2))
                except:
                    pass
            if s.get("omega") and s.get("L"):
                try:
                    g_vals.append(s["omega"]**2 * s["L"])
                except:
                    pass
        if g_vals:
            g_mean = sum(g_vals)/len(g_vals)
            for s in states:
                if s.get("g") is None:
                    s["g"] = g_mean

    # Infer omega from two (x,v) pairs if available
    pairs = [s for s in states if s.get("x") is not None and s.get("v") is not None]
    if len(pairs) >= 2:
        # Use first two pairs
        x1, v1 = pairs[0]["x"], pairs[0]["v"]
        x2, v2 = pairs[1]["x"], pairs[1]["v"]
        denom = (x2**2 - x1**2)
        if abs(denom) > 1e-12:
            omega_sq = (v1**2 - v2**2) / denom
            if omega_sq > 0:
                omega_est = math.sqrt(omega_sq)
                for s in states:
                    if s.get("omega") is None:
                        s["omega"] = omega_est
                # if pendulum, propagate g via L if present
                if system_type == 2:
                    for s in states:
                        if s.get("L") and s.get("g") is None:
                            s["g"] = omega_est**2 * s["L"]

def finalize_state(s, system_type):
    # iterative passes to fill derived quantities
    for _ in range(8):
        # omega from T or f or k/m or g/L
        if s.get("omega") is None:
            if s.get("T") is not None:
                try: s["omega"] = 2*math.pi / s["T"]
                except: pass
            elif s.get("f") is not None:
                s["omega"] = 2*math.pi * s["f"]
            elif system_type == 1 and s.get("k") is not None and s.get("m") is not None and s["m"] != 0:
                s["omega"] = math.sqrt(s["k"]/s["m"])
            elif system_type == 2 and s.get("L") is not None and s.get("g") is not None and s["L"] != 0:
                s["omega"] = math.sqrt(s["g"]/s["L"])

        # default sigmaF
        if s.get("sigmaF") is None:
            if system_type == 1:
                if s.get("k") is not None and s.get("x") is not None:
                    s["sigmaF"] = s["k"] * s["x"]
                elif s.get("m") is not None and s.get("g") is not None:
                    s["sigmaF"] = s["m"] * s["g"]
            else:
                if s.get("m") is not None and s.get("g") is not None:
                    s["sigmaF"] = s["m"] * s["g"]

        # k from sigmaF and x (spring)
        if system_type == 1 and s.get("k") is None and s.get("sigmaF") is not None and s.get("x") is not None:
            if abs(s["x"]) > 1e-12:
                s["k"] = s["sigmaF"] / s["x"]

        # k/m <-> omega
        if system_type == 1:
            if s.get("k") is None and s.get("m") and s.get("omega"):
                s["k"] = s["m"] * s["omega"]**2
            if s.get("m") is None and s.get("k") and s.get("omega") and s["omega"] != 0:
                s["m"] = s["k"] / s["omega"]**2

        # A from E or Vmax or v,x
        if s.get("A") is None:
            if s.get("v") is not None and s.get("x") is not None and s.get("omega") is not None and s["omega"] != 0:
                val = (s["v"]**2)/(s["omega"]**2) + s["x"]**2
                if val >= 0:
                    s["A"] = math.sqrt(val)
            if s.get("A") is None and s.get("Vmax") is not None and s.get("omega") is not None and s["omega"] != 0:
                s["A"] = s["Vmax"] / s["omega"]
            if s.get("A") is None and s.get("a_max") is not None and s.get("omega") is not None and s["omega"] != 0:
                s["A"] = s["a_max"] / (s["omega"]**2)
            if s.get("A") is None and s.get("E") is not None and s.get("k") is not None and s["k"]>0:
                s["A"] = math.sqrt(2*s["E"]/s["k"])
            if s.get("A") is None and s.get("E") is not None and s.get("m") is not None and s.get("omega") is not None and s["omega"]!=0:
                s["A"] = math.sqrt(2*s["E"]/(s["m"]*s["omega"]**2))

        # v if missing from v formula
        if s.get("v") is None:
            if s.get("omega") is not None and s.get("A") is not None and s.get("x") is not None:
                val = s["omega"]**2 * (s["A"]**2 - s["x"]**2)
                if val >= 0:
                    s["v"] = math.sqrt(val)
            elif s.get("KE") is not None and s.get("m") is not None and s["m"]>0:
                s["v"] = math.sqrt(2*s["KE"]/s["m"])
            elif s.get("E") is not None and s.get("PE") is not None and s.get("m") is not None and s["m"]>0:
                s["v"] = math.sqrt(max(0, 2*(s["E"] - s["PE"])/s["m"]))

        # Vmax, a_max
        if s.get("Vmax") is None and s.get("omega") is not None and s.get("A") is not None:
            s["Vmax"] = s["omega"]*s["A"]
        if s.get("a_max") is None and s.get("omega") is not None and s.get("A") is not None:
            s["a_max"] = s["omega"]**2 * s["A"]

        # a instantaneous
        if s.get("a") is None and s.get("omega") is not None and s.get("x") is not None:
            s["a"] = - s["omega"]**2 * s["x"]

        # energies
        if s.get("PE") is None and s.get("k") is not None and s.get("x") is not None:
            s["PE"] = 0.5 * s["k"] * s["x"]**2
        if s.get("E") is None and s.get("k") is not None and s.get("A") is not None:
            s["E"] = 0.5 * s["k"] * s["A"]**2
        if s.get("KE") is None and s.get("m") is not None and s.get("v") is not None:
            s["KE"] = 0.5 * s["m"] * s["v"]**2
        if s.get("KE") is None and s.get("E") is not None and s.get("PE") is not None:
            s["KE"] = s["E"] - s["PE"]
        if s.get("E") is None and s.get("KE") is not None and s.get("PE") is not None:
            s["E"] = s["KE"] + s["PE"]

        # f,T
        if s.get("f") is None and s.get("omega") is not None:
            s["f"] = s["omega"] / (2*math.pi)
        if s.get("T") is None and s.get("f") is not None and s["f"] != 0:
            s["T"] = 1.0 / s["f"]

    return s

def process_states(raw_states, system="spring"):
    # normalize raw inputs: convert numbers/"="/"None"
    system_type = 1 if str(system).lower().startswith("s") else 2
    states = []
    global_mem = {}

    # first pass: parse tokens and remember numeric values
    for raw in raw_states:
        s = {}
        for k in NUMERIC_KEYS:
            v = raw.get(k, None)
            parsed = to_num(v)
            if parsed == "=":
                # placeholder: if global_mem has value, use it; else keep "=" for now
                s[k] = global_mem.get(k, "=")
            else:
                s[k] = parsed
                if isinstance(parsed, (int,float)):
                    global_mem[k] = parsed
        states.append(s)

    # backfill: if any remaining "=" placeholders and later states define them, find later definition
    known = {k:v for k,v in global_mem.items()}
    for k in NUMERIC_KEYS:
        if known.get(k) is None:
            for s in states:
                val = s.get(k)
                if val is not None and val != "=":
                    known[k] = val
                    break
    # replace "=" in states with known if possible
    for s in states:
        for k in NUMERIC_KEYS:
            if s.get(k) == "=":
                s[k] = known.get(k)

    # infer parameters using helper
    infer_global_params(states, system_type)

    # finalize each state
    results = []
    for s in states:
        s_copy = copy.deepcopy(s)
        if system_type == 2 and s_copy.get("g") is None:
            s_copy["g"] = 9.81
        s_fin = finalize_state(s_copy, system_type)
        # format floats or None
        out = {}
        for k in NUMERIC_KEYS:
            out[k] = fmt(s_fin.get(k))
        out["warnings"] = warn_if_inconsistent(s_fin, system_type)
        out["explanation"] = explain_physics(system_type)
        results.append(out)
    return results

# quick CLI test when run directly
if __name__ == "__main__":
    import json
    print("Run small CLI test")
    payload = {
        "system":"spring",
        "states":[
            {"m":2,"x":0.10},
            {"m":0.5,"x":"="}
        ]
    }
    res = process_states(payload["states"], payload["system"])
    print(json.dumps(res, indent=2))
