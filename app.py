# app.py
from flask import Flask, render_template, request, jsonify
from shm_solver import process_states
import json

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/solve", methods=["POST"])
def solve():
    payload = request.get_json(force=True)
    system = payload.get("system", "spring")
    raw_states = payload.get("states", [])
    # Allow keys with Thai ΣF or 'sigmaF'
    normalized_states = []
    for st in raw_states:
        norm = {}
        for k in st:
            # normalize possible ΣF key variants
            if k in ("ΣF","SigmaF","sigmaF"):
                norm["sigmaF"] = st[k]
            else:
                norm[k] = st[k]
        normalized_states.append(norm)
    results = process_states(normalized_states, system)
    return jsonify({"system":system, "results": results})

if __name__ == "__main__":
    app.run(debug=True)
