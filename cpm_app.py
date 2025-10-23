import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from collections import defaultdict

st.set_page_config(page_title="Critical Path Method Analysis", layout="wide")

# --- Formatter angka: 0 -> "0", 1.5 -> "1,5", 1.50000 -> "1,5" ---
def fmt(x):
    try:
        xv = float(x)
        if xv.is_integer():
            return str(int(xv))
        return f"{xv:g}".replace(".", ",")
    except Exception:
        return str(x)

# --------------------- UI STYLES ---------------------
st.markdown("""
<style>
    .main-header {
        text-align: center;
        color: #1f77b4;
        padding: 20px;
        background: linear-gradient(90deg, #e3f2fd 0%, #bbdefb 100%);
        border-radius: 10px;
        margin-bottom: 30px;
    }
    .metric-card {
        background: white;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
    }
    .tutorial-box {
        background: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        border-left: 4px solid #1f77b4;
        margin: 10px 0;
    }
    .formula-table {
        font-size: 0.85em;
        margin: 10px 0;
    }
    .example-box {
        background: #fff3cd;
        padding: 12px;
        border-radius: 6px;
        border-left: 4px solid #ffc107;
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header"><h1>📊 Critical Path Method (CPM) Analysis By Athila</h1><p>Network Diagram & Schedule Management</p></div>', unsafe_allow_html=True)

# --------------------- DATA ---------------------
if 'activities' not in st.session_state:
    st.session_state.activities = [
        {"Activity": "A", "Initial Node": 1, "Final Node": 2, "Duration": 2},
        {"Activity": "B", "Initial Node": 2, "Final Node": 3, "Duration": 2},
        {"Activity": "C", "Initial Node": 2, "Final Node": 4, "Duration": 3},
        {"Activity": "D", "Initial Node": 2, "Final Node": 5, "Duration": 4},
        {"Activity": "E", "Initial Node": 3, "Final Node": 6, "Duration": 2},
        {"Activity": "F", "Initial Node": 4, "Final Node": 6, "Duration": 3},
        {"Activity": "G", "Initial Node": 5, "Final Node": 7, "Duration": 6},
        {"Activity": "H", "Initial Node": 6, "Final Node": 8, "Duration": 2},
        {"Activity": "I", "Initial Node": 6, "Final Node": 7, "Duration": 5},
        {"Activity": "J", "Initial Node": 7, "Final Node": 8, "Duration": 1},
        {"Activity": "K", "Initial Node": 8, "Final Node": 9, "Duration": 2},
    ]


with st.sidebar:
    st.header("📝 Input Data")
    st.info("Edit tabel di bawah untuk mengubah data aktivitas")
    df = pd.DataFrame(st.session_state.activities)
    edited_df = st.data_editor(
        df, num_rows="dynamic", use_container_width=True,
        column_config={
            "Activity": st.column_config.TextColumn("Activity", required=True),
            "Initial Node": st.column_config.NumberColumn("Initial Node", required=True, min_value=1, format="%g"),
            "Final Node": st.column_config.NumberColumn("Final Node", required=True, min_value=1, format="%g"),
            "Duration": st.column_config.NumberColumn("Duration", required=True, min_value=0, format="%g"),
        }
    )
    if st.button("🔄 Update & Calculate", type="primary", use_container_width=True):
        st.session_state.activities = edited_df.to_dict('records')
        st.rerun()

    st.markdown("---")
    
    # --------------------- TUTORIAL PDM ---------------------
    st.markdown("### 📖 Tutorial PDM")
    
    with st.expander("🔹 1. Rumus Dasar CPM", expanded=False):
        st.markdown("""
        <div class="tutorial-box">
        <table class="formula-table" style="width:100%; border-collapse: collapse;">
        <thead>
            <tr style="background-color: #e3f2fd;">
                <th style="padding: 8px; border: 1px solid #ddd;">Simbol</th>
                <th style="padding: 8px; border: 1px solid #ddd;">Arti</th>
                <th style="padding: 8px; border: 1px solid #ddd;">Rumus</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>ES</strong><br/>(Early Start)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">Waktu paling awal kegiatan bisa dimulai</td>
                <td style="padding: 8px; border: 1px solid #ddd;">= <strong>EF</strong> kegiatan sebelumnya <strong>terbesar</strong></td>
            </tr>
            <tr style="background-color: #f8f9fa;">
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>EF</strong><br/>(Early Finish)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">Waktu paling awal kegiatan bisa selesai</td>
                <td style="padding: 8px; border: 1px solid #ddd;">= <strong>ES + Duration</strong></td>
            </tr>
            <tr>
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>LF</strong><br/>(Late Finish)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">Waktu paling akhir kegiatan bisa selesai tanpa tunda proyek</td>
                <td style="padding: 8px; border: 1px solid #ddd;">= <strong>LS</strong> kegiatan sesudahnya <strong>terkecil</strong></td>
            </tr>
            <tr style="background-color: #f8f9fa;">
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>LS</strong><br/>(Late Start)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">Waktu paling akhir kegiatan bisa dimulai tanpa tunda proyek</td>
                <td style="padding: 8px; border: 1px solid #ddd;">= <strong>LF - Duration</strong></td>
            </tr>
            <tr>
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>TS</strong><br/>(Total Slack)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">Waktu longgar total (total float)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">= <strong>LS - ES</strong><br/>atau <strong>LF - EF</strong></td>
            </tr>
            <tr style="background-color: #f8f9fa;">
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>FF</strong><br/>(Free Float)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">Waktu longgar bebas (tanpa pengaruh ke successor)</td>
                <td style="padding: 8px; border: 1px solid #ddd;">= <strong>ES(successor) - EF(current)</strong></td>
            </tr>
        </tbody>
        </table>
        </div>
        """, unsafe_allow_html=True)
    
    with st.expander("🔹 2. Aturan Gabungan Node (Multiple Predecessors)", expanded=False):
        st.markdown("""
        <div class="tutorial-box">
        <p><strong>Aturan:</strong> Kalau sebuah aktivitas punya <strong>dua atau lebih pendahulu</strong>, maka ES-nya diambil dari <strong>EF terbesar</strong> dari semua pendahulunya.</p>
        </div>
        
        <div class="example-box">
        <p><strong>📌 Contoh:</strong></p>
        
        <table style="width:100%; border-collapse: collapse; font-size: 0.9em; margin: 10px 0;">
        <thead>
            <tr style="background-color: #fff3cd;">
                <th style="padding: 6px; border: 1px solid #ddd;">Aktivitas</th>
                <th style="padding: 6px; border: 1px solid #ddd;">Duration</th>
                <th style="padding: 6px; border: 1px solid #ddd;">Predecessor</th>
            </tr>
        </thead>
        <tbody>
            <tr><td style="padding: 6px; border: 1px solid #ddd;">A</td><td style="padding: 6px; border: 1px solid #ddd;">3</td><td style="padding: 6px; border: 1px solid #ddd;">–</td></tr>
            <tr><td style="padding: 6px; border: 1px solid #ddd;">B</td><td style="padding: 6px; border: 1px solid #ddd;">4</td><td style="padding: 6px; border: 1px solid #ddd;">A</td></tr>
            <tr><td style="padding: 6px; border: 1px solid #ddd;">C</td><td style="padding: 6px; border: 1px solid #ddd;">6</td><td style="padding: 6px; border: 1px solid #ddd;">A</td></tr>
            <tr><td style="padding: 6px; border: 1px solid #ddd;">D</td><td style="padding: 6px; border: 1px solid #ddd;">5</td><td style="padding: 6px; border: 1px solid #ddd;">B, C</td></tr>
        </tbody>
        </table>
        
        <p><strong>Perhitungan:</strong></p>
        <ul style="font-size: 0.9em; line-height: 1.6;">
        <li><strong>A:</strong> ES = 0, EF = 0 + 3 = <strong>3</strong></li>
        <li><strong>B:</strong> ES = 3, EF = 3 + 4 = <strong>7</strong></li>
        <li><strong>C:</strong> ES = 3, EF = 3 + 6 = <strong>9</strong></li>
        <li><strong>D:</strong> Karena D punya dua pendahulu (B dan C), ambil <strong>EF terbesar</strong>:<br/>
        ES(D) = max(EF<sub>B</sub>, EF<sub>C</sub>) = max(7, 9) = <strong>9</strong><br/>
        EF(D) = 9 + 5 = <strong>14</strong></li>
        </ul>
        
        <p style="margin-top: 10px; padding: 8px; background-color: #d1ecf1; border-radius: 4px; font-size: 0.85em;">
        💡 <strong>Kesimpulan:</strong> Aktivitas D harus menunggu sampai aktivitas C selesai (yang lebih lama) sebelum bisa dimulai.
        </p>
        </div>
        """, unsafe_allow_html=True)
    
    with st.expander("🔹 3. Aturan Percabangan Node (Multiple Successors)", expanded=False):
        st.markdown("""
        <div class="tutorial-box">
        <p><strong>Aturan:</strong> Kalau sebuah aktivitas punya <strong>dua atau lebih penerus</strong>, maka LF-nya diambil dari <strong>LS terkecil</strong> dari semua penerusnya.</p>
        </div>
        
        <div class="example-box">
        <p><strong>📌 Contoh (Lanjutan dari contoh sebelumnya):</strong></p>
        
        <p>Misal dari aktivitas D ada dua cabang:</p>
        
        <table style="width:100%; border-collapse: collapse; font-size: 0.9em; margin: 10px 0;">
        <thead>
            <tr style="background-color: #fff3cd;">
                <th style="padding: 6px; border: 1px solid #ddd;">Aktivitas</th>
                <th style="padding: 6px; border: 1px solid #ddd;">Duration</th>
                <th style="padding: 6px; border: 1px solid #ddd;">Predecessor</th>
            </tr>
        </thead>
        <tbody>
            <tr><td style="padding: 6px; border: 1px solid #ddd;">E</td><td style="padding: 6px; border: 1px solid #ddd;">4</td><td style="padding: 6px; border: 1px solid #ddd;">D</td></tr>
            <tr><td style="padding: 6px; border: 1px solid #ddd;">F</td><td style="padding: 6px; border: 1px solid #ddd;">6</td><td style="padding: 6px; border: 1px solid #ddd;">D</td></tr>
        </tbody>
        </table>
        
        <p><strong>Backward Pass (dari akhir proyek):</strong></p>
        <p>Misalkan LF proyek = 20</p>
        
        <ul style="font-size: 0.9em; line-height: 1.6;">
        <li><strong>E:</strong> LF = 20, LS = 20 - 4 = <strong>16</strong></li>
        <li><strong>F:</strong> LF = 20, LS = 20 - 6 = <strong>14</strong></li>
        <li><strong>D:</strong> Karena D punya dua successor (E dan F), ambil <strong>LS terkecil</strong>:<br/>
        LF(D) = min(LS<sub>E</sub>, LS<sub>F</sub>) = min(16, 14) = <strong>14</strong><br/>
        LS(D) = 14 - 5 = <strong>9</strong></li>
        </ul>
        
        <p style="margin-top: 10px; padding: 8px; background-color: #d1ecf1; border-radius: 4px; font-size: 0.85em;">
        💡 <strong>Kesimpulan:</strong> Aktivitas D harus selesai maksimal di waktu 14 agar aktivitas F (yang paling ketat) tidak terlambat.
        </p>
        </div>
        """, unsafe_allow_html=True)
    
    with st.expander("🔹 4. Interpretasi Slack & Critical Path", expanded=False):
        st.markdown("""
        <div class="tutorial-box">
        <table style="width:100%; border-collapse: collapse; font-size: 0.9em;">
        <thead>
            <tr style="background-color: #e3f2fd;">
                <th style="padding: 8px; border: 1px solid #ddd;">Nilai TS</th>
                <th style="padding: 8px; border: 1px solid #ddd;">Arti</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>0</strong></td>
                <td style="padding: 8px; border: 1px solid #ddd;">Aktivitas <strong>kritis</strong> (harus tepat waktu, tidak ada kelonggaran)</td>
            </tr>
            <tr style="background-color: #f8f9fa;">
                <td style="padding: 8px; border: 1px solid #ddd;"><strong>&gt; 0</strong></td>
                <td style="padding: 8px; border: 1px solid #ddd;">Aktivitas <strong>non-kritis</strong> (masih ada kelonggaran waktu)</td>
            </tr>
        </tbody>
        </table>
        </div>
        
        <div class="example-box">
        <p><strong>🎯 Critical Path (Jalur Kritis):</strong></p>
        <ul style="font-size: 0.9em; line-height: 1.6;">
        <li>Jalur yang terdiri dari aktivitas-aktivitas dengan <strong>TS = 0</strong></li>
        <li>Jalur yang menentukan <strong>durasi minimum proyek</strong></li>
        <li>Jika ada keterlambatan di jalur ini, <strong>seluruh proyek akan terlambat</strong></li>
        <li>Ditampilkan dengan <strong>warna merah</strong> pada diagram</li>
        </ul>
        
        <p style="margin-top: 10px;"><strong>📊 Total Slack (TS) vs Free Float (FF):</strong></p>
        <ul style="font-size: 0.9em; line-height: 1.6;">
        <li><strong>Total Slack (TS):</strong> Berapa lama aktivitas bisa ditunda tanpa menunda proyek</li>
        <li><strong>Free Float (FF):</strong> Berapa lama aktivitas bisa ditunda tanpa mempengaruhi aktivitas berikutnya</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with st.expander("🔹 5. Contoh Visualisasi PDM", expanded=False):
        st.markdown("""
        <div class="tutorial-box">
        <p><strong>Format Box Aktivitas di Diagram PDM:</strong></p>
        </div>
        
        <div style="text-align: center; margin: 15px 0;">
        <svg width="200" height="140" xmlns="http://www.w3.org/2000/svg">
            <!-- Box -->
            <rect x="10" y="10" width="180" height="120" rx="8" fill="#cce5ff" stroke="black" stroke-width="2"/>
            
            <!-- Activity Name -->
            <text x="100" y="50" font-size="20" font-weight="bold" text-anchor="middle" fill="black">A</text>
            <text x="100" y="70" font-size="12" text-anchor="middle" fill="black" font-style="italic">D=3</text>
            
            <!-- ES (top-left) -->
            <text x="30" y="30" font-size="10" font-weight="bold" text-anchor="middle" fill="blue">ES</text>
            <text x="30" y="42" font-size="10" font-weight="bold" text-anchor="middle" fill="blue">0</text>
            
            <!-- EF (top-right) -->
            <text x="170" y="30" font-size="10" font-weight="bold" text-anchor="middle" fill="blue">EF</text>
            <text x="170" y="42" font-size="10" font-weight="bold" text-anchor="middle" fill="blue">3</text>
            
            <!-- LS (bottom-left) -->
            <text x="30" y="105" font-size="10" font-weight="bold" text-anchor="middle" fill="green">LS</text>
            <text x="30" y="117" font-size="10" font-weight="bold" text-anchor="middle" fill="green">0</text>
            
            <!-- LF (bottom-right) -->
            <text x="170" y="105" font-size="10" font-weight="bold" text-anchor="middle" fill="green">LF</text>
            <text x="170" y="117" font-size="10" font-weight="bold" text-anchor="middle" fill="green">3</text>
            
            <!-- FF label (above box) -->
            <rect x="75" y="-5" width="50" height="18" rx="4" fill="#c8e6c9" stroke="black" stroke-width="1"/>
            <text x="100" y="8" font-size="9" text-anchor="middle" fill="black">FF: 0</text>
            
            <!-- TS label (below box) -->
            <rect x="75" y="132" width="50" height="18" rx="4" fill="#fff9c4" stroke="black" stroke-width="1"/>
            <text x="100" y="145" font-size="9" text-anchor="middle" fill="black">TS: 0</text>
        </svg>
        </div>
        
        <div class="example-box" style="font-size: 0.85em;">
        <p><strong>Keterangan:</strong></p>
        <ul style="line-height: 1.6;">
        <li><strong style="color: blue;">ES (Early Start)</strong> - Pojok kiri atas</li>
        <li><strong style="color: blue;">EF (Early Finish)</strong> - Pojok kanan atas</li>
        <li><strong style="color: green;">LS (Late Start)</strong> - Pojok kiri bawah</li>
        <li><strong style="color: green;">LF (Late Finish)</strong> - Pojok kanan bawah</li>
        <li><strong>D (Duration)</strong> - Di tengah, di bawah nama aktivitas</li>
        <li><strong>FF (Free Float)</strong> - Label di atas box (background hijau)</li>
        <li><strong>TS (Total Slack)</strong> - Label di bawah box (background kuning)</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("---")
    st.markdown("### 📖 Legend")
    st.markdown("**ES** = Early Start  •  **EF** = Early Finish")
    st.markdown("**LS** = Late Start   •  **LF** = Late Finish")
    st.markdown("**TS** = Total Slack  •  **FF** = Free Float")

# --------------------- CPM CORE ---------------------
def calculate_cpm(activities):
    G = nx.DiGraph()
    for act in activities:
        G.add_edge(act['Initial Node'], act['Final Node'],
                   activity=act['Activity'], duration=act['Duration'])

    nodes = sorted(G.nodes())
    start_node, end_node = min(nodes), max(nodes)

    # Forward pass (event times)
    # ----- Topological order untuk perhitungan event time -----
    topo = list(nx.topological_sort(G))

    # ----- Forward pass: Earliest Event Time (TE) per node -----
    TE = {n: 0 for n in topo}
    for u in topo:
        for v in G.successors(u):
            d = G[u][v]['duration']
            TE[v] = max(TE[v], TE[u] + d)

    # Project duration = earliest time di node akhir (bukan EF[end_node])
    project_duration = TE[end_node]

    # ----- Backward pass: Latest Event Time (TL) per node -----
    TL = {n: project_duration for n in topo}
    for u in reversed(topo):
        if G.out_degree(u) == 0:
            TL[u] = project_duration
        else:
            TL[u] = min(TL[v] - G[u][v]['duration'] for v in G.successors(u))

    # ----- Activity metrics (pakai TE/TL yang bener) -----
    metrics = {}
    for a in activities:
        start, end, d = a['Initial Node'], a['Final Node'], a['Duration']
        es = TE[start]
        ef = es + d
        lf = TL[end]
        ls = lf - d
        ts = ls - es
        ff = TE[end] - ef
        metrics[a['Activity']] = dict(
            ES=es, EF=ef, LS=ls, LF=lf, TS=ts, FF=ff,
            start_node=start, end_node=end, duration=d
        )
    return metrics, project_duration, G

def find_all_paths(G, start, end, path=None):
    if path is None: path = []
    path = path + [start]
    if start == end: return [path]
    res = []
    for n in G.successors(start):
        if n not in path:
            res.extend(find_all_paths(G, n, end, path))
    return res

def get_path_info(G, path):
    acts, dur = [], 0
    for i in range(len(path)-1):
        e = G[path[i]][path[i+1]]
        acts.append(e['activity']); dur += e['duration']
    return acts, dur

# --------------------- AOA LAYOUT (lane-based) ---------------------
def hierarchical_layout_aoa(G):
    topo = list(nx.topological_sort(G))
    TE = {n: 0 for n in topo}
    for u in topo:
        for v in G.successors(u):
            TE[v] = max(TE[v], TE[u] + G[u][v].get('duration', 0))

    uniq_te = sorted(set(TE.values()))
    col_of = {t: i for i, t in enumerate(uniq_te)}
    node_col = {n: col_of[TE[n]] for n in topo}

    columns = defaultdict(list)
    for n in topo:
        columns[node_col[n]].append(n)

    x_spacing, y_spacing = 3.5, 2.0
    pos = {}

    # siblings map
    siblings = defaultdict(list)
    for p in G.nodes():
        succs = list(G.successors(p))
        if len(succs) >= 2:
            siblings[p] = succs[:]

    # initial order per column
    for c in columns:
        def bary(n):
            sc = [node_col[s] for s in G.successors(n)]
            return (sum(sc) / len(sc)) if sc else c
        columns[c].sort(key=lambda n: (bary(n), n))
    for c in sorted(columns):
        for j, n in enumerate(columns[c]):
            pos[n] = (c * x_spacing, 0.0)  # Y diset oleh lane

    # lanes
    parent_lane = {}
    start_node = min(topo)
    for n in topo:
        if n == start_node:
            parent_lane[n] = 0
        elif n not in parent_lane:
            preds = list(G.predecessors(n))
            parent_lane[n] = parent_lane.get(preds[0], 0) if preds else 0

    # sebar anak multi-branch: -1,+1,-2,+2,...
    for n in topo:
        succs = list(G.successors(n))
        if len(succs) >= 2:
            def child_priority(s):
                hits_sib = sum(1 for t in G.successors(s) if t in succs)
                downstream = len(nx.descendants(G, s))
                return (-hits_sib, -downstream, s)
            succs_sorted = sorted(succs, key=child_priority)
            offs, d = [], 1
            while len(offs) < len(succs_sorted):
                offs.append(-d)
                if len(offs) < len(succs_sorted):
                    offs.append(+d)
                d += 1
            for i, s in enumerate(succs_sorted):
                parent_lane[s] = parent_lane[n] + offs[i]

    # ahead bias tipis untuk node yang terima/lempar ke sibling
    for c in sorted(columns):
        by_par = defaultdict(list)
        for n in columns[c]:
            key = tuple(sorted(G.predecessors(n)))
            by_par[key].append(n)
        for sibs in by_par.values():
            if len(sibs) <= 1:
                continue
            for n in sibs:
                x, y = pos[n]
                succ = list(G.successors(n))
                in_from_sib = any(p in sibs for p in G.predecessors(n))
                out_to_sib = any(s in sibs for s in succ)
                ax = 0.0
                if in_from_sib:
                    ax += 0.40 * x_spacing
                if out_to_sib:
                    ax += 0.15 * x_spacing
                if any(G.out_degree(s) == 0 for s in succ):
                    ax += 0.12 * x_spacing
                pos[n] = (x + ax, y)

    # set Y dari lane + jitter halus biar gak numpuk
    for c in sorted(columns):
        groups = defaultdict(list)
        for n in columns[c]:
            groups[tuple(sorted(G.predecessors(n)))].append(n)
        for sibs in groups.values():
            if len(sibs) == 1:
                n = sibs[0]
                x, _ = pos[n]
                pos[n] = (x, -parent_lane.get(n, 0) * y_spacing)
                continue
            def vscore(n):
                succ = list(G.successors(n))
                out_to_sib = sum(1 for s in succ if s in sibs)
                in_from_sib = sum(1 for p in G.predecessors(n) if p in sibs)
                branch = len(succ)
                return (2 * out_to_sib) - (2 * in_from_sib) + 0.5 * branch
            sibs.sort(key=lambda n: (-vscore(n), n))
            k = len(sibs)
            for i, n in enumerate(sibs):
                x, _ = pos[n]
                base_y = -parent_lane.get(n, 0) * y_spacing
                jitter = (i - (k - 1) / 2) * 0.25
                pos[n] = (x, base_y + jitter)

    # dorong turun/naik kalau ada “long jump” (gap kolom besar)
    for n in G.nodes():
        x, y = pos[n]
        succ = list(G.successors(n))
        pred = list(G.predecessors(n))
        out_jumps = sum(max(0, (node_col[s] - node_col[n])) for s in succ)
        in_jumps = sum(max(0, (node_col[n] - node_col[p])) for p in pred)
        y += 0.40 * (out_jumps - in_jumps)
        pos[n] = (x, y)

    # ekstra drop untuk anak dari long jump (mis. 2→6)
    for n in G.nodes():
        x, y = pos[n]
        for p in G.predecessors(n):
            col_gap = node_col[n] - node_col[p]
            if col_gap >= 2:
                y -= 0.8 * y_spacing * (col_gap - 1)
        pos[n] = (x, y)

    # align sumber paling kiri & muara paling kanan
    first_col, last_col = min(columns), max(columns)
    for n in G.nodes():
        x, y = pos[n]
        if G.in_degree(n) == 0:
            pos[n] = (first_col * x_spacing, y)
        if G.out_degree(n) == 0:
            pos[n] = (last_col * x_spacing, y)

    return pos, parent_lane


# --------------------- MAIN LOGIC ---------------------
try:
    activity_metrics, project_duration, G = calculate_cpm(st.session_state.activities)
    critical_activities = [a for a, m in activity_metrics.items() if m['TS'] == 0]

    nodes = sorted(G.nodes())
    start_node, end_node = min(nodes), max(nodes)
    all_paths = find_all_paths(G, start_node, end_node)
    path_info = []
    for p in all_paths:
        acts, dur = get_path_info(G, p)
        path_info.append({
            'Path': ' → '.join(acts),
            'Duration': dur,
            'Critical': 'Yes' if dur == project_duration else 'No'
        })

    # --------------------- TABS ---------------------
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary", "🗺️ AOA Network Diagram", "📐 PDM Diagram", "📋 Detailed Metrics"])

    # --------------------- SUMMARY ---------------------
    with tab1:
        st.markdown("### 🎯 Project Summary")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("🕒 Project Duration", f"{project_duration} days")
        with col2:
            st.metric("📌 Critical Activities", len(critical_activities))
        with col3:
            st.metric("🔗 Total Activities", len(st.session_state.activities))

        st.markdown("#### 🔴 Critical Path Activities")
        st.info(' → '.join(critical_activities))

        st.markdown("#### 📍 All Possible Paths")
        path_df = pd.DataFrame(path_info)
        def highlight_critical_path(row):
            return ['background-color: #ffcccc']*len(row) if row['Critical']=='Yes' else ['']*len(row)
        st.dataframe(
            path_df.style
                .apply(highlight_critical_path, axis=1)
                .format(fmt),
            use_container_width=True,
            hide_index=True
        )
    # --------------------- AOA DIAGRAM ---------------------
    with tab2:
        st.markdown("### 🗺️ Activity-on-Arrow (AOA) Network Diagram")
        fig, ax = plt.subplots(figsize=(16, 10))
        pos, lanes = hierarchical_layout_aoa(G)  # lanes kebaca, gak wajib dipakai

        # Identify critical edges
        critical_edges = set()
        for i in range(len(critical_activities) - 1):
            curr = critical_activities[i]
            nxt = critical_activities[i + 1]
            for a in st.session_state.activities:
                if a['Activity'] == curr:
                    start_node_curr = a['Final Node']
                if a['Activity'] == nxt:
                    start_node_nxt = a['Initial Node']
            if start_node_curr == start_node_nxt:
                # Find edge
                for a in st.session_state.activities:
                    if a['Activity'] == curr:
                        u, v = a['Initial Node'], a['Final Node']
                        critical_edges.add((u, v))
                for a in st.session_state.activities:
                    if a['Activity'] == nxt:
                        u, v = a['Initial Node'], a['Final Node']
                        critical_edges.add((u, v))

        # Draw all edges with labels
        for u, v, data in G.edges(data=True):
            x1, y1 = pos[u]; x2, y2 = pos[v]
            is_critical = (u, v) in critical_edges
            color = 'red' if is_critical else 'gray'
            lw = 3 if is_critical else 1.5
            ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, zorder=1)
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            label = f"{data['activity']}\n({data['duration']})"
            ax.text(mx, my + 0.3, label, fontsize=11, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', edgecolor='black', linewidth=1.2),
                    fontweight='bold' if is_critical else 'normal', zorder=3)

        # Draw nodes
        for n in G.nodes():
            x, y = pos[n]
            color = '#90EE90' if n == start_node else '#FFB6C6' if n == end_node else 'lightblue'
            circle = plt.Circle((x, y), 0.5, color=color, ec='black', linewidth=2.5, zorder=2)
            ax.add_patch(circle)
            ax.text(x, y, str(n), fontsize=16, fontweight='bold', ha='center', va='center', zorder=4)

        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='red', lw=3, label='Critical Path'),
            Line2D([0], [0], color='gray', lw=1.5, label='Non-Critical Path'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#90EE90', markersize=12,
                   markeredgecolor='black', markeredgewidth=2, label='Start Node', linestyle='None'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#FFB6C6', markersize=12,
                   markeredgecolor='black', markeredgewidth=2, label='End Node', linestyle='None'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.95)
        ax.set_title("Activity-on-Arrow (AOA) Network Diagram", fontsize=18, fontweight='bold', pad=20)
        ax.axis('equal'); ax.axis('off'); plt.tight_layout()
        st.pyplot(fig)

    # --------------------- PDM DIAGRAM ---------------------
    with tab3:
        st.markdown("### 📐 Precedence Diagramming Method (PDM)")

        # lane assignment by activity metrics
        def act_lane(act_name):
            m = activity_metrics[act_name]
            if m['TS'] == 0: return 0  # critical lane
            # non-critical: group by ES
            es_val = m['ES']
            return 1 + (es_val % 3)

        # build ES-level groups
        es_levels = defaultdict(list)
        for a in st.session_state.activities:
            m = activity_metrics[a['Activity']]
            es_levels[m['ES']].append(a['Activity'])

        sorted_es = sorted(es_levels.keys())
        x_spacing, y_spacing = 4.5, 2.5
        pdm_pos = {}

        # position per ES level
        for level_idx, es in enumerate(sorted_es):
            # kelompokkan per lane dalam level yang sama
            by_lane = defaultdict(list)
            for a in es_levels[es]:
                by_lane[act_lane(a)].append(a)
            for ln, arr in by_lane.items():
                arr.sort()  # biar deterministik
                k = len(arr)
                for j, a in enumerate(arr):
                    # spread kecil sekitar pusat lane, supaya gak nabrak
                    offset = (j - (k - 1) / 2) * 1.1
                    pdm_pos[a] = (level_idx * x_spacing, -ln * y_spacing + offset)

        # START/FINISH positions follow critical lane
        first_es = min(sorted_es)
        last_es  = max(sorted_es)

        level_first_ys = [pdm_pos[a][1] for a in es_levels[first_es]]
        level_last_ys  = [pdm_pos[a][1] for a in es_levels[last_es]]

        def median(vals): 
            s = sorted(vals); n=len(s)
            return s[n//2] if n%2==1 else 0.5*(s[n//2-1]+s[n//2])

        start_y = median(level_first_ys)
        finish_x = max(x for x,_ in pdm_pos.values()) + x_spacing
        finish_y = median(level_last_ys)

        fig, ax = plt.subplots(figsize=(18, 12))

        # START box
        start_x = min(x for x,_ in pdm_pos.values()) - x_spacing
        start_box = FancyBboxPatch((start_x-0.8, start_y-0.7), 1.6, 1.4,
                                   boxstyle="round,pad=0.1", edgecolor='black', facecolor='#90EE90', linewidth=3)
        ax.add_patch(start_box)
        ax.text(start_x, start_y, "START", ha='center', va='center', fontsize=15, fontweight='bold')

        # FINISH box
        finish_box = FancyBboxPatch((finish_x-0.8, finish_y-0.7), 1.6, 1.4,
                                    boxstyle="round,pad=0.1", edgecolor='black', facecolor='#FFB6C6', linewidth=3)
        ax.add_patch(finish_box)
        ax.text(finish_x, finish_y+0.15, "FINISH", ha='center', va='center', fontsize=15, fontweight='bold')
        ax.text(finish_x, finish_y-0.35, f"{project_duration} days", ha='center', va='center', fontsize=12, style='italic')

        # draw activity boxes
        from matplotlib.lines import Line2D
        for a in st.session_state.activities:
            name=a['Activity']; m=activity_metrics[name]
            crit = name in critical_activities
            x,y = pdm_pos[name]
            BOX_W = 2.4   # lebar box aktivitas (default sebelumnya 1.8)
            BOX_H = 1.8   # tinggi box aktivitas (default sebelumnya 1.3)
            CORNER_PAD = 0.10  # paddings sudut untuk FancyBboxPatch
            FONT_MAIN = 16     # nama aktivitas
            FONT_DUR = 11      # "D=..."
            FONT_CORNER = 10   # ES/EF/LS/LF
            FONT_TAG = 9       # TS/FF

            # Offset teks di dalam box (relatif ke BOX_W/BOX_H)
            # kiri/kanan atas/bawah (pojok)
            OFF_X = BOX_W/2 - 0.2
            OFF_Y = BOX_H/2 - 0.17

            # label TS/FF di luar box
            TAG_GAP = 0.30  # jarak label TS/FF dari sisi box

            # Panah masuk/keluar box
            ARROW_INSET = BOX_W/2   # jarak ujung panah dari pusat menuju tepi box

            # Perbesar jarak antar level/lane kalau kotaknya makin gede
            box = FancyBboxPatch((x-BOX_W/2, y-BOX_H/2), BOX_W, BOX_H,
                     boxstyle=f"round,pad={CORNER_PAD}",
                     edgecolor='black',
                     facecolor=('#ffcccc' if crit else '#cce5ff'),
                     linewidth=(3 if crit else 1.8))
            ax.add_patch(box)

            # Nama & Durasi (posisi relatif supaya proporsional saat diperbesar)
            ax.text(x, y+0.08, name, ha='center', va='center',
                    fontsize=FONT_MAIN, fontweight='bold')
            ax.text(x, y-0.25, f"D={m['duration']}", ha='center', va='center',
                    fontsize=FONT_DUR, style='italic')

            # ES/EF/LS/LF sudut (pakai OFF_X/OFF_Y)
            ax.text(x-OFF_X, y+OFF_Y, f"ES\n{fmt(m['ES'])}", ha='center', va='center',
                    fontsize=FONT_CORNER, color='blue', fontweight='bold')
            ax.text(x+OFF_X, y+OFF_Y, f"EF\n{fmt(m['EF'])}", ha='center', va='center',
                    fontsize=FONT_CORNER, color='blue', fontweight='bold')
            ax.text(x-OFF_X, y-OFF_Y, f"LS\n{fmt(m['LS'])}", ha='center', va='center',
                    fontsize=FONT_CORNER, color='green', fontweight='bold')
            ax.text(x+OFF_X, y-OFF_Y, f"LF\n{fmt(m['LF'])}", ha='center', va='center',
                    fontsize=FONT_CORNER, color='green', fontweight='bold')

            # TS bawah & FF atas (pakai TAG_GAP dari sisi box)
            ax.text(x, y-(BOX_H/2 + TAG_GAP), f"TS: {fmt(m['TS'])}", ha='center', va='center',
                    fontsize=FONT_TAG, bbox=dict(boxstyle='round,pad=0.35', facecolor='#fff9c4', edgecolor='black', linewidth=1.2))
            ax.text(x, y+(BOX_H/2 + TAG_GAP), f"FF: {fmt(m['FF'])}", ha='center', va='center',
                    fontsize=FONT_TAG, bbox=dict(boxstyle='round,pad=0.35', facecolor='#c8e6c9', edgecolor='black', linewidth=1.2))
                        


        # arrows between activities
        drawn_arrows=set()
        for a in st.session_state.activities:
            act=a['Activity']; end_n=a['Final Node']
            for nxt in st.session_state.activities:
                if nxt['Initial Node'] != end_n: continue
                nxt_act = nxt['Activity']
                if (act, nxt_act) in drawn_arrows: continue
                if act not in pdm_pos or nxt_act not in pdm_pos: continue
                # skip start/finish edges here (handled below)
                if act=="START" or nxt_act=="FINISH": continue

                x1,y1=pdm_pos[act]; x2,y2=pdm_pos[nxt_act]
                is_crit = (act in critical_activities and nxt_act in critical_activities)
                color = 'red' if is_crit else 'gray'
                lw = 3 if is_crit else 1.5
                dy=abs(y2-y1); rad = 0.25 if dy>2.5 else 0.15
                arr = FancyArrowPatch((x1 + ARROW_INSET, y1), (x2 - ARROW_INSET, y2),
                      arrowstyle='->', mutation_scale=28,
                      color=color, linewidth=lw, connectionstyle=f"arc3,rad={rad}")
                ax.add_patch(arr); drawn_arrows.add((act,nxt_act))

        # arrows from START to first-level activities
        first_node = min(a['Initial Node'] for a in st.session_state.activities)
        for a in st.session_state.activities:
            if a['Initial Node'] == first_node:
                act=a['Activity']
                if act in pdm_pos:
                    x,y = pdm_pos[act]
                    arr = FancyArrowPatch((start_x+0.8, start_y),(x-0.9,y),
                                          arrowstyle='->', mutation_scale=25,
                                          color='gray', linewidth=2, connectionstyle="arc3,rad=0.15")
                    ax.add_patch(arr)

        # arrows from last-level activities to FINISH
        last_node = max(a['Final Node'] for a in st.session_state.activities)
        for a in st.session_state.activities:
            if a['Final Node'] == last_node:
                act=a['Activity']
                if act in pdm_pos:
                    x,y = pdm_pos[act]
                    is_crit = act in critical_activities
                    arr = FancyArrowPatch((x+0.9,y),(finish_x-0.8, finish_y),
                                          arrowstyle='->', mutation_scale=25,
                                          color=('red' if is_crit else 'gray'),
                                          linewidth=(3 if is_crit else 2),
                                          connectionstyle="arc3,rad=0.15")
                    ax.add_patch(arr)

        # legend & title
        from matplotlib.patches import Rectangle
        legend_elements = [
            Line2D([0],[0],color='red',lw=3,label='Critical Path'),
            Line2D([0],[0],color='gray',lw=1.5,label='Non-Critical Path'),
            Rectangle((0,0),1,1,fc='#ffcccc',ec='black',lw=2,label='Critical Activity'),
            Rectangle((0,0),1,1,fc='#cce5ff',ec='black',lw=1.5,label='Non-Critical Activity'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.95)

        info_text=("ES = Early Start (Top Left) | EF = Early Finish (Top Right)\n"
                   "LS = Late Start (Bottom Left) | LF = Late Finish (Bottom Right)\n"
                   "TS = Total Slack (Below) | FF = Free Float (Above) | D = Duration")
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10, va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.95, pad=0.8))

        ax.set_title("Precedence Diagramming Method (PDM) with Complete Time Analysis", fontsize=18, fontweight='bold', pad=20)
        ax.axis('equal'); ax.axis('off'); plt.tight_layout()
        st.pyplot(fig)

    # --------------------- METRICS ---------------------
    with tab4:
        st.markdown("### 📋 Complete Activity Metrics")
        rows=[]
        for a in st.session_state.activities:
            name=a['Activity']; m=activity_metrics[name]; crit=(name in critical_activities)
            rows.append({'Activity':name,'Start Node':m['start_node'],'End Node':m['end_node'],
                         'Duration':m['duration'],'ES':m['ES'],'EF':m['EF'],'LS':m['LS'],'LF':m['LF'],
                         'Total Slack':m['TS'],'Free Float':m['FF'],'Critical':'✓' if crit else ''})
        metrics_df = pd.DataFrame(rows)
        def highlight_critical(row):
            return ['background-color: #ffcccc']*len(row) if row['Critical']=='✓' else ['']*len(row)
        st.dataframe(
            metrics_df.style
                .apply(highlight_critical, axis=1)
                .format(fmt, subset=["Duration","ES","EF","LS","LF","Total Slack","Free Float"]),
            use_container_width=True,
            hide_index=True
        )
        st.markdown("---")
        c1,c2 = st.columns(2)
        with c1:
            st.download_button("📥 Download Metrics (CSV)", data=metrics_df.to_csv(index=False),
                               file_name="cpm_metrics.csv", mime="text/csv", use_container_width=True)
        with c2:
            path_df = pd.DataFrame(path_info)
            st.download_button("📥 Download Path Analysis (CSV)", data=path_df.to_csv(index=False),
                               file_name="path_analysis.csv", mime="text/csv", use_container_width=True)

except Exception as e:
    st.error(f"❌ Error in calculation: {str(e)}")
    st.info("Please check your input data. Make sure all nodes are properly connected and form a valid network.")
    st.exception(e)

st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
  <p><strong>📚 Critical Path Method (CPM) Analysis Tool</strong></p>
  <p>Based on Project Time Management & Network Scheduling Principles</p>
  <p style='font-size: 0.9em;'>ES=Early Start | EF=Early Finish | LS=Late Start | LF=Late Finish | TS=Total Slack | FF=Free Float</p>
</div>
""", unsafe_allow_html=True)
