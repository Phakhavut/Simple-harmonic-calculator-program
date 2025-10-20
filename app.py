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
<!DOCTYPE html>
<html lang="th">
<head>
<meta charset="utf-8"/>
<title>Simple-harmonic-calculator-program </title>
<script src="https://cdn.plot.ly/plotly-2.30.0.min.js"></script>
<style>
  body{font-family:Arial,Helvetica,sans-serif;margin:18px;background:#f7f9fb;color:#111}
  h1{margin-bottom:6px}
  .box{background:white;padding:12px;border-radius:8px;box-shadow:0 1px 4px rgba(0,0,0,0.08);margin-bottom:14px}
  h3{margin-top:0}
  label{display:inline-block;width:140px;font-weight:bold;font-size:0.9em}
  input[type="text"]{width:70px;padding:4px;border:1px solid #ddd;border-radius:4px;margin-right:4px}
  button{padding:8px 12px;margin:4px;background:#1677ff;color:white;border:none;border-radius:6px;cursor:pointer}
  button:hover{background:#0d5eff}
  .state-row{margin-bottom:12px;padding:10px;border:1px solid #e6eef6;border-radius:6px;background:#f9fcff}
  .state-row h4{margin:0 0 8px 0;font-size:1.1em}
  .vars-container{display:flex;flex-wrap:wrap}
  .var-group{display:inline-block;margin:4px 6px;vertical-align:top;min-width:160px}
  .warning{color:#b91c1c;font-weight:700}
  table{width:100%;border-collapse:collapse;margin-top:8px}
  td,th{border:1px solid #e6eef6;padding:6px;text-align:left}
  th{background:#f0f8ff}
  #result, #plots{display:none}
  #plots .plot-div{height:360px;margin-top:12px}
  select{width:120px;padding:4px}
  .hidden{display:none !important}
</style>
</head>
<body>
  <h1>Simple-harmonic-calculator-program</h1>

  <div class="box">
    <h3>🔧 คลิกเพื่อทดลองตัวอย่าง</h3>
    <button onclick="loadExample('spring1')">สปริง: m=2kg → m=0.5kg (x คง)</button>
    <button onclick="loadExample('spring2')">สปริง: m=0.05kg, x_eq=0.02m, A=0.025m</button>
    <button onclick="loadExample('pendulum1')">ลูกตุ้ม: L=2m, T=2.5s → L=8m</button>
    <button onclick="clearAll()">ล้างทั้งหมด</button>
  </div>

  <div class="box">
    <h3>📝 ใส่ข้อมูล </h3>
    <div>
      <label>ระบบ:</label>
      <select id="system" onchange="updateVisibility()">
        <option value="spring">สปริง (Spring)</option>
        <option value="pendulum">ลูกตุ้ม (Pendulum)</option>
      </select>
    </div>
    <div style="margin:12px 0">
      <label>จำนวนสภาวะ:</label>
      <select id="numStates" onchange="generateStates()">
        <option value="1">1</option>
        <option value="2" selected>2</option>
        <option value="3">3</option>
        <option value="4">4</option>
      </select>
      <button onclick="addState()">เพิ่มสภาวะ</button>
      <button onclick="removeLastState()">ลบสภาวะสุดท้าย</button>
    </div>
    <div id="statesContainer"></div>
    <br>
    <button onclick="calculate()">คำนวณ</button>
    <button onclick="exportJson()">Export เป็น JSON (ดูโค้ด)</button>
  </div>

  <div id="result" class="box"></div>

<script>
let currentNumStates = 2;
const varLabels = {
  m: 'มวล m (kg)', k: 'สปริง k (N/m)', L: 'ความยาว L (m)', g: 'แรงโน้มถ่วง g (m/s²)',
  sigmaF: 'ΣF (N)', x: 'การกระจัด x (m)', v: 'ความเร็ว v (m/s)', Vmax: 'V_max (m/s)',
  a_max: 'a_max (m/s²)', a: 'ความเร่ง a (m/s²)', KE: 'KE (J)', PE: 'PE (J)',
  E: 'พลังงานรวม E (J)', A: 'แอมพลิจูด A (m)', omega: 'ω (rad/s)',
  f: 'ความถี่ f (Hz)', T: 'รอบ T (s)'
};
const allVars = ['m','k','L','g','sigmaF','x','v','Vmax','a_max','a','KE','PE','E','A','omega','f','T'];

function generateStates() {
  const num = parseInt(document.getElementById('numStates').value);
  currentNumStates = num;
  const container = document.getElementById('statesContainer');
  container.innerHTML = '';
  for (let i = 0; i < num; i++) {
    const row = createStateRow(i + 1);
    container.appendChild(row);
  }
  updateVisibility();
}

function createStateRow(index) {
  const div = document.createElement('div');
  div.className = 'state-row';
  div.innerHTML = `<h4>สภาวะ ${index}</h4><div class="vars-container"></div>`;
  const varsContainer = div.querySelector('.vars-container');
  allVars.forEach(key => {
    const group = document.createElement('div');
    group.className = 'var-group';
    const label = document.createElement('label');
    label.textContent = varLabels[key] || key;
    const input = document.createElement('input');
    input.type = 'text';  // เปลี่ยนเป็น text เพื่อให้กรอก = ได้
    input.id = `state${index}_${key}`;
    group.appendChild(label);
    group.appendChild(input);
    varsContainer.appendChild(group);
  });
  return div;
}

function updateVisibility() {
  const system = document.getElementById('system').value;
  const isSpring = system === 'spring';
  allVars.forEach(key => {
    for (let i = 1; i <= currentNumStates; i++) {
      const input = document.getElementById(`state${i}_${key}`);
      const group = input ? input.parentElement : null;
      if (group) {
        if ((key === 'k' && !isSpring) || (key === 'L' && isSpring)) {
          group.classList.add('hidden');
        } else {
          group.classList.remove('hidden');
        }
      }
    }
  });
}

function addState() {
  currentNumStates++;
  const container = document.getElementById('statesContainer');
  container.appendChild(createStateRow(currentNumStates));
  updateVisibility();
}

function removeLastState() {
  if (currentNumStates <= 1) return;
  const container = document.getElementById('statesContainer');
  container.removeChild(container.lastChild);
  currentNumStates--;
}

function getInputValues() {
  const system = document.getElementById('system').value;
  const states = [];
  for (let i = 1; i <= currentNumStates; i++) {
    const state = {};
    allVars.forEach(key => {
      const input = document.getElementById(`state${i}_${key}`);
      if (input && input.value !== '') {
        state[key] = input.value.trim();
      }
    });
    states.push(state);
  }
  return {system, states};
}

function calculate() {
  const payload = getInputValues();
  fetch('/solve', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(payload)
  })
  .then(res => res.json())
  .then(data => {
    const resultDiv = document.getElementById('result');
    if (data.error) {
      resultDiv.innerHTML = `<p class="warning">เกิดข้อผิดพลาด: ${data.error}</p>`;
    } else {
      let html = `<h3>ผลลัพธ์ (${data.system})</h3>`;
      data.results.forEach((s, idx) => {
        html += `<h4>สภาวะ ${idx+1}</h4><table><tr>`;
        for (let key of allVars) {
          html += `<th>${varLabels[key] || key}</th>`;
        }
        html += `</tr><tr>`;
        for (let key of allVars) {
          html += `<td>${s[key] !== null ? s[key] : ''}</td>`;
        }
        html += `</tr></table>`;
        if (s.warnings && s.warnings.length > 0) {
          html += `<p class="warning">คำเตือน: ${s.warnings.join('; ')}</p>`;
        }
        html += `<pre>${s.explanation}</pre>`;
      });
      resultDiv.innerHTML = html;
      resultDiv.style.display = 'block';
    }
  })
  .catch(err => {
    document.getElementById('result').innerHTML = `<p class="warning">เกิดข้อผิดพลาด: ${err}</p>`;
  });
}

function clearAll() {
  document.getElementById('statesContainer').innerHTML = '';
  generateStates();
}

function loadExample(name) {
  let system, values;
  if (name === 'spring1') {
    system = 'spring';
    values = [{m:2,x:0.1,g:9.81},{m:0.5,x:'='}];
  } else if (name === 'spring2') {
    system = 'spring';
    values = [{m:0.05,x:0.02,g:9.81},{m:'=',A:0.025,x:0}];
  } else {
    system = 'pendulum';
    values = [{L:2,T:2.5,g:9.81},{L:8,T:'='}];
  }
  document.getElementById('system').value = system;
  document.getElementById('numStates').value = values.length;
  currentNumStates = values.length;
  generateStates();
  values.forEach((s, idx) => {
    for (let k in s) {
      const input = document.getElementById(`state${idx+1}_${k}`);
      if (input) input.value = s[k];
    }
  });
  updateVisibility();
}

function exportJson() {
  const payload = getInputValues();
  alert(JSON.stringify(payload, null, 2));
}

// initialize
generateStates();
</script>
</body>
</html>


