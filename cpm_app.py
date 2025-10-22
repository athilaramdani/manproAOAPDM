import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from collections import defaultdict

st.set_page_config(page_title="Critical Path Method Analysis", layout="wide")

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
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header"><h1>📊 Critical Path Method (CPM) Analysis</h1><p>Network Diagram & Schedule Management</p></div>', unsafe_allow_html=True)

# --------------------- DATA ---------------------
if 'activities' not in st.session_state:
    st.session_state.activities = [
        {"Activity": "A", "Initial Node": 1, "Final Node": 2, "Duration": 3},
        {"Activity": "B", "Initial Node": 2, "Final Node": 3, "Duration": 4},
        {"Activity": "C", "Initial Node": 2, "Final Node": 6, "Duration": 6},
        {"Activity": "D", "Initial Node": 3, "Final Node": 4, "Duration": 9},
        {"Activity": "E", "Initial Node": 3, "Final Node": 5, "Duration": 3},
        {"Activity": "F", "Initial Node": 4, "Final Node": 5, "Duration": 6},
        {"Activity": "G", "Initial Node": 5, "Final Node": 6, "Duration": 9},
        {"Activity": "H", "Initial Node": 4, "Final Node": 7, "Duration": 5},
        {"Activity": "I", "Initial Node": 5, "Final Node": 7, "Duration": 8},
        {"Activity": "J", "Initial Node": 6, "Final Node": 7, "Duration": 2},
    ]

with st.sidebar:
    st.header("📝 Input Data")
    st.info("Edit tabel di bawah untuk mengubah data aktivitas")
    df = pd.DataFrame(st.session_state.activities)
    edited_df = st.data_editor(
        df, num_rows="dynamic", use_container_width=True,
        column_config={
            "Activity": st.column_config.TextColumn("Activity", required=True),
            "Initial Node": st.column_config.NumberColumn("Initial Node", required=True, min_value=1),
            "Final Node": st.column_config.NumberColumn("Final Node", required=True, min_value=1),
            "Duration": st.column_config.NumberColumn("Duration", required=True, min_value=0),
        }
    )
    if st.button("🔄 Update & Calculate", type="primary", use_container_width=True):
        st.session_state.activities = edited_df.to_dict('records')
        st.rerun()

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
    col_of = {t:i for i,t in enumerate(uniq_te)}
    node_col = {n: col_of[TE[n]] for n in topo}

    columns = defaultdict(list)
    for n in topo: columns[node_col[n]].append(n)

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
            return (sum(sc)/len(sc)) if sc else c
        columns[c].sort(key=lambda n: (bary(n), n))
    for c in sorted(columns):
        col_nodes = columns[c]
        for j, n in enumerate(col_nodes):
            pos[n] = (c * x_spacing, 0.0)  # Y akan di-set oleh lane

    # lanes
    parent_lane = {}
    start_node = min(topo)
    for n in topo:
        if n == start_node:
            parent_lane[n] = 0
        elif n not in parent_lane:
            preds = list(G.predecessors(n))
            parent_lane[n] = parent_lane.get(preds[0], 0) if preds else 0

    # spread multi-children: -1, +1, -2, +2, ...
    for n in topo:
        succs = list(G.successors(n))
        if len(succs) >= 2:
            def child_priority(s):
                hits_sib = sum(1 for t in G.successors(s) if t in succs)
                downstream = len(nx.descendants(G, s))
                return (-hits_sib, -downstream, s)
            succs_sorted = sorted(succs, key=child_priority)
            offs = []
            d = 1
            while len(offs) < len(succs_sorted):
                offs.append(-d)
                if len(offs) < len(succs_sorted): offs.append(+d)
                d += 1
            for i, s in enumerate(succs_sorted):
                parent_lane[s] = parent_lane[n] + offs[i]

    # ahead bias (sedikit) untuk yang menerima dari sibling
    for c in sorted(columns):
        by_par = defaultdict(list)
        for n in columns[c]:
            key = tuple(sorted(G.predecessors(n)))
            by_par[key].append(n)
        for sibs in by_par.values():
            if len(sibs) <= 1: continue
            for n in sibs:
                x, y = pos[n]
                succ = list(G.successors(n))
                in_from_sib = any(p in sibs for p in G.predecessors(n))
                out_to_sib = any(s in sibs for s in succ)
                ax = 0.0
                if in_from_sib: ax += 0.40 * x_spacing
                if out_to_sib:  ax += 0.15 * x_spacing
                if any(G.out_degree(s) == 0 for s in succ): ax += 0.12 * x_spacing
                pos[n] = (x + ax, y)

    # Y from lane (murni) + jitter tipis
    for c in sorted(columns):
        col_nodes = columns[c]
        groups = defaultdict(list)
        for n in col_nodes:
            groups[tuple(sorted(G.predecessors(n)))] .append(n)
        for sibs in groups.values():
            if len(sibs) == 1:
                n = sibs[0]
                x,_ = pos[n]
                pos[n] = (x, -parent_lane.get(n,0)*y_spacing)
                continue
            def vscore(n):
                succ = list(G.successors(n))
                out_to_sib = sum(1 for s in succ if s in sibs)
                in_from_sib = sum(1 for p in G.predecessors(n) if p in sibs)
                branch = len(succ)
                return (2*out_to_sib) - (2*in_from_sib) + 0.5*branch
            sibs.sort(key=lambda n:(-vscore(n), n))
            k = len(sibs)
            for i,n in enumerate(sibs):
                x,_ = pos[n]
                base_y = -parent_lane.get(n,0)*y_spacing
                jitter = (i - (k-1)/2) * 0.25
                pos[n] = (x, base_y + jitter)

    # push long-jump children a bit lower/higher
    for n in G.nodes():
        x,y = pos[n]
        succ = list(G.successors(n)); pred = list(G.predecessors(n))
        out_jumps = sum(max(0,(node_col[s]-node_col[n])) for s in succ)
        in_jumps  = sum(max(0,(node_col[n]-node_col[p])) for p in pred)
        y += 0.40*(out_jumps - in_jumps)
        pos[n] = (x,y)

    # extra drop for nodes that came from long jump (e.g., 2->6)
    for n in G.nodes():
        x,y = pos[n]
        for p in G.predecessors(n):
            col_gap = node_col[n]-node_col[p]
            if col_gap >= 2:
                y -= 0.8*y_spacing*(col_gap-1)
        pos[n] = (x,y)

    # force left-most sources, right-most sinks aligned to extreme columns
    first_col, last_col = min(columns), max(columns)
    for n in G.nodes():
        x,y = pos[n]
        if G.in_degree(n)==0:  pos[n]=(first_col*x_spacing, y)
        if G.out_degree(n)==0: pos[n]=(last_col*x_spacing,  y)

    return pos, parent_lane

# --------------------- MAIN CALC ---------------------
try:
    activity_metrics, project_duration, G = calculate_cpm(st.session_state.activities)
    nodes = sorted(G.nodes())
    start_node, end_node = min(nodes), max(nodes)

    all_paths = find_all_paths(G, start_node, end_node)
    path_info = []
    for p in all_paths:
        acts, dur = get_path_info(G, p)
        path_info.append({'path': ' → '.join(map(str,p)),
                          'activities': ' → '.join(acts),
                          'duration': dur})
    path_info.sort(key=lambda x: x['duration'], reverse=True)
    critical_path_info = path_info[0]
    critical_activities = critical_path_info['activities'].split(' → ')

    tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary", "🗺️ AOA Diagram", "📦 PDM Diagram", "📋 Detailed Metrics"])

    # --------------------- SUMMARY ---------------------
    with tab1:
        c1,c2,c3,c4 = st.columns(4)
        with c1: st.markdown('<div class="metric-card">', unsafe_allow_html=True); st.metric("🎯 Project Duration", f"{project_duration} days"); st.markdown('</div>', unsafe_allow_html=True)
        with c2: st.markdown('<div class="metric-card">', unsafe_allow_html=True); st.metric("📌 Total Activities", len(st.session_state.activities)); st.markdown('</div>', unsafe_allow_html=True)
        with c3: st.markdown('<div class="metric-card">', unsafe_allow_html=True); st.metric("🔴 Critical Activities", len(critical_activities)); st.markdown('</div>', unsafe_allow_html=True)
        with c4: st.markdown('<div class="metric-card">', unsafe_allow_html=True); st.metric("🛤️ Total Paths", len(path_info)); st.markdown('</div>', unsafe_allow_html=True)

        st.markdown("---")
        L,R = st.columns([1,1])
        with L:
            st.markdown("### 🎯 Critical Path")
            st.success(f"**Activities:** {critical_path_info['activities']}")
            st.info(f"**Nodes:** {critical_path_info['path']}")
            st.warning(f"**Duration:** {critical_path_info['duration']} days")
            crit_df = pd.DataFrame([{'Activity':a,'Duration':activity_metrics[a]['duration'],'Total Slack':activity_metrics[a]['TS']} for a in critical_activities])
            st.dataframe(crit_df, use_container_width=True, hide_index=True)
        with R:
            st.markdown("### 🛤️ All Possible Paths")
            for i,p in enumerate(path_info,1):
                is_crit = i==1
                icon = "🔴" if is_crit else "🔵"
                with st.expander(f"{icon} Path {i} - {p['duration']} days {'(CRITICAL)' if is_crit else ''}"):
                    st.write(f"**Nodes:** {p['path']}")
                    st.write(f"**Activities:** {p['activities']}")
                    st.write(f"**Duration:** {p['duration']} days")

    # --------------------- AOA ---------------------
    with tab2:
        st.markdown("### 🗺️ Activity-on-Arrow (AOA) Network Diagram")
        fig, ax = plt.subplots(figsize=(20,12))
        pos, lanes = hierarchical_layout_aoa(G)

        # nodes
        for n in G.nodes():
            x,y = pos[n]
            circ = plt.Circle((x,y), 0.35, color='lightblue', ec='black', linewidth=2.5, zorder=4)
            ax.add_patch(circ)
            ax.text(x, y, str(n), ha='center', va='center', fontsize=16, fontweight='bold', zorder=5)

        # edges + labels (mid-curve)
        edge_count = defaultdict(int)
        for act in st.session_state.activities:
            u,v = act['Initial Node'], act['Final Node']
            a = act['Activity']; d = act['Duration']
            if u not in pos or v not in pos: continue
            crit = a in critical_activities
            color = 'red' if crit else 'gray'
            lw = 3.5 if crit else 2

            x1,y1 = pos[u]; x2,y2 = pos[v]
            dx,dy = x2-x1, y2-y1
            L = np.hypot(dx,dy); 
            if L < 1e-6: continue
            nx_, ny_ = dx/L, dy/L
            sx, sy = x1+nx_*0.35, y1+ny_*0.35
            ex, ey = x2-nx_*0.35, y2-ny_*0.35

            edge_key=(u,v); edge_count[edge_key]+=1; k=edge_count[edge_key]
            x_spacing_local=3.5
            col_gap=max(1,int(round(abs(dx)/x_spacing_local)))
            rad=(0.18+0.08*(col_gap-1))
            if k%2==0: rad=-rad
            if abs(dy)<0.3: rad*=1.2

            arrow = FancyArrowPatch((sx,sy),(ex,ey), arrowstyle='->', mutation_scale=25,
                                    color=color, linewidth=lw, connectionstyle=f"arc3,rad={rad}", zorder=2)
            ax.add_patch(arrow)

            mx,my=(x1+x2)/2,(y1+y2)/2
            k_off=0.5
            ox,oy= -ny_*rad*k_off, nx_*rad*k_off
            lx,ly= mx+ox, my+oy
            bbox_color = '#ffcccc' if crit else '#ffffcc'
            ax.text(lx, ly, f"{a} = {d}", fontsize=11, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor=bbox_color, edgecolor='black', linewidth=1.5),
                    ha='center', zorder=6)

        ax.set_title("Activity-on-Arrow (AOA) Network Diagram\n(Semua aktivitas dari tabel ditampilkan dengan panah melengkung)", fontsize=18, fontweight='bold', pad=20)
        ax.axis('equal'); ax.axis('off')
        from matplotlib.lines import Line2D
        ax.legend([Line2D([0],[0],color='red',lw=3.5,label='Critical Path'),
                   Line2D([0],[0],color='gray',lw=2,label='Non-Critical Path')],
                   ['Critical Path','Non-Critical Path'], loc='upper right', fontsize=12)
        plt.tight_layout()
        st.pyplot(fig)

        with st.expander("🔍 Verifikasi: Semua Aktivitas yang Digambar"):
            verify_df = pd.DataFrame([{
                'Activity': a['Activity'],
                'From Node': a['Initial Node'],
                'To Node':   a['Final Node'],
                'Duration':  a['Duration'],
                'Status':    '🔴 Critical' if a['Activity'] in critical_activities else '⚪ Normal'
            } for a in st.session_state.activities])
            st.dataframe(verify_df, use_container_width=True, hide_index=True)

    # --------------------- PDM ---------------------
    with tab3:
        st.markdown("### 📦 Precedence Diagramming Method (PDM) - Activity on Node")
        fig, ax = plt.subplots(figsize=(22,14))

        # group by ES level
        es_levels = defaultdict(list)
        for a in st.session_state.activities:
            es = activity_metrics[a['Activity']]['ES']
            es_levels[es].append(a['Activity'])
        sorted_es = sorted(es_levels.keys())

        # map lane from AOA by initial node
        def act_lane(act_name):
            row = next(r for r in st.session_state.activities if r["Activity"]==act_name)
            return lanes.get(row["Initial Node"], 0)

        x_spacing, y_spacing = 4.5, 3.0
        pdm_pos = {}
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
            box = FancyBboxPatch((x-0.9, y-0.65), 1.8, 1.3, boxstyle="round,pad=0.08",
                                 edgecolor='black', facecolor=('#ffcccc' if crit else '#cce5ff'),
                                 linewidth=(3 if crit else 1.8))
            ax.add_patch(box)
            ax.text(x, y+0.08, name, ha='center', va='center', fontsize=16, fontweight='bold')
            ax.text(x, y-0.25, f"D={m['duration']}", ha='center', va='center', fontsize=11, style='italic')

            cs=9
            ax.text(x-0.7, y+0.48, f"ES\n{m['ES']}", ha='center', va='center', fontsize=cs, color='blue', fontweight='bold')
            ax.text(x+0.7, y+0.48, f"EF\n{m['EF']}", ha='center', va='center', fontsize=cs, color='blue', fontweight='bold')
            ax.text(x-0.7, y-0.48, f"LS\n{m['LS']}", ha='center', va='center', fontsize=cs, color='green', fontweight='bold')
            ax.text(x+0.7, y-0.48, f"LF\n{m['LF']}", ha='center', va='center', fontsize=cs, color='green', fontweight='bold')

            ax.text(x, y-0.95, f"TS: {m['TS']}", ha='center', va='center',
                    fontsize=9, bbox=dict(boxstyle='round,pad=0.35', facecolor='#fff9c4', edgecolor='black', linewidth=1.2))
            ax.text(x, y+0.95, f"FF: {m['FF']}", ha='center', va='center',
                    fontsize=9, bbox=dict(boxstyle='round,pad=0.35', facecolor='#c8e6c9', edgecolor='black', linewidth=1.2))

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
                arr = FancyArrowPatch((x1+0.9,y1),(x2-0.9,y2), arrowstyle='->', mutation_scale=25,
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
        st.dataframe(metrics_df.style.apply(highlight_critical, axis=1), use_container_width=True, hide_index=True)

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
