import math

def parse_float(val, prev=None):
    """แปลงค่าเป็น float ถ้าว่างหรือ '=' ให้ใช้ค่า prev"""
    if val.strip() == "" or val.strip() == "=":
        return prev
    try:
        return float(val)
    except:
        return None

def calculate_states(system, states):
    """
    system: "spring" หรือ "pendulum"
    states: list ของ dict ของแต่ละ state
    คืนค่า list ของ dict ที่มีค่าเต็ม
    """
    results = []
    prev = {}
    g_default = 9.81

    for s in states:
        state = {}
        # อ่านค่าพื้นฐาน
        state['m'] = parse_float(s.get('m'), prev.get('m'))
        state['k'] = parse_float(s.get('k'), prev.get('k'))
        state['L'] = parse_float(s.get('L'), prev.get('L'))
        state['g'] = parse_float(s.get('g'), prev.get('g', g_default))
        state['sigmaF'] = parse_float(s.get('sigmaF'), prev.get('sigmaF'))

        # ค่าตำแหน่งและความเร็ว
        state['x'] = parse_float(s.get('x'), prev.get('x'))
        state['v'] = parse_float(s.get('v'), prev.get('v'))

        # คำนวณค่าที่เหลือ
        # สำหรับสปริง
        if system == "spring":
            if state['x'] is not None and state['k'] is not None and state['m'] is not None:
                state['Vmax'] = state['v'] if state['v'] is not None else math.sqrt(state['k']/state['m'] * state['x']**2)
                state['a_max'] = math.sqrt(state['k']/state['m']) * abs(state['x'])
                state['a'] = state['a_max'] * (-1 if state['x']>0 else 1)
                state['KE'] = 0.5 * state['m'] * state['v']**2 if state['v'] is not None else 0.5 * state['m'] * (state['Vmax']**2)
                state['PE'] = 0.5 * state['k'] * state['x']**2
                state['E'] = state['KE'] + state['PE']
                state['A'] = abs(state['x'])
                state['omega'] = math.sqrt(state['k']/state['m'])
                state['f'] = state['omega']/(2*math.pi)
                state['T'] = 1/state['f'] if state['f'] != 0 else None

        # สำหรับลูกตุ้ม
        elif system == "pendulum":
            if state['L'] is not None:
                state['omega'] = math.sqrt(state['g']/state['L'])
                state['f'] = state['omega']/(2*math.pi)
                state['T'] = 1/state['f'] if state['f'] != 0 else None
                if state['x'] is not None:
                    state['Vmax'] = state['v'] if state['v'] is not None else state['omega']*abs(state['x'])
                    state['a_max'] = state['omega']**2 * abs(state['x'])
                    state['a'] = -state['omega']**2 * state['x']
                    state['KE'] = 0.5 * state['m'] * state['v']**2 if state['v'] is not None else 0.5 * state['m'] * (state['Vmax']**2)
                    state['PE'] = 0.5 * state['m'] * state['g'] * state['x']**2
                    state['E'] = state['KE'] + state['PE']
                    state['A'] = abs(state['x'])

        # กรณีอื่นๆ ใส่ None
        for k in ['Vmax','a_max','a','KE','PE','E','A','omega','f','T']:
            if k not in state:
                state[k] = None

        prev = state
        results.append(state)
    return results
