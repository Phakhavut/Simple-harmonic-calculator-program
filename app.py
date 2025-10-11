# app.py
from flask import Flask, render_template, request, jsonify
from shm_solver import process_states
import json
import traceback  # สำหรับ debug

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/solve", methods=["POST"])
def solve():
    try:
        payload = request.get_json(force=True)
        if not payload:
            raise ValueError("Invalid or empty JSON payload")
        
        system = payload.get("system", "spring")
        raw_states = payload.get("states", [])
        if not isinstance(raw_states, list) or len(raw_states) == 0:
            raise ValueError("No 'states' array in payload")
        
        # Normalize keys (handle ΣF variants, including Thai ΣF)
        normalized_states = []
        for st in raw_states:
            if not isinstance(st, dict):
                raise ValueError("Each state must be a dict")
            norm = {}
            for k in st:
                if k in ("ΣF", "SigmaF", "sigmaF", "ΣF"):  # รองรับตัวอักษรไทย
                    norm["sigmaF"] = st[k]
                else:
                    norm[k] = st[k]
            normalized_states.append(norm)
        
        results = process_states(normalized_states, system)
        return jsonify({"system": system, "results": results})
    
    except Exception as e:
        error_msg = str(e)
        # Log full traceback to console for debug
        print("ERROR in /solve:", error_msg)
        print(traceback.format_exc())
        return jsonify({"error": error_msg, "results": []}), 400

if __name__ == "__main__":
    app.run(debug=True)
