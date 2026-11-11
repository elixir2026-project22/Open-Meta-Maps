"""Module detail page: sharded script data (works from file://) plus a static viewer page."""
import json
import ui
from core import INPUTS

_RAW={}

def stub_far(g,m):
    """Segments reaching a metabolite far outside the module (source-map long-range links) become short dashed stubs, so one stray endpoint cannot stretch the figure."""
    import math
    segs=g['s']
    if not segs:return g
    L=[math.hypot(s[0]-s[2],s[1]-s[3]) for s in segs];med=sorted(L)[len(L)//2];thr=max(2200,10*med)
    xs=sorted(l[1] for l in g['l']);ys=sorted(l[2] for l in g['l'])
    cx,cy=(xs[len(xs)//2],ys[len(ys)//2]) if xs else (segs[0][0],segs[0][1])
    mets={x[0]:x for x in g['m']};used=set();ghosts=[];g['far']=0
    by={(x[2],x[3]):x[0] for x in g['m']}
    for s,l in zip(segs,L):
        if l<=thr:
            for p in((s[0],s[1]),(s[2],s[3])):
                if p in by:used.add(by[p])
            continue
        d0=math.hypot(s[0]-cx,s[1]-cy);d1=math.hypot(s[2]-cx,s[3]-cy);far_end=1 if d1>=d0 else 0
        nx,ny,fx,fy=(s[0],s[1],s[2],s[3]) if far_end else (s[2],s[3],s[0],s[1])
        stub=min(l,max(4*med,150));ux,uy=(fx-nx)/l,(fy-ny)/l;sx,sy=round(nx+ux*stub,1),round(ny+uy*stub,1)
        met=mets.get(by.get((fx,fy)))
        if far_end:s[2],s[3]=sx,sy
        else:s[0],s[1]=sx,sy
        s[4:8]=[None]*4;s.append(1)
        ghosts.append([(by.get((fx,fy)) or 'x')+'~',met[1] if met else '',sx,sy,0,round(sx+8,1),round(sy-8,1)]);g['far']+=1
    for s in segs:
        if len(s)<12:s.append(0)
    g['m']=[x for x in g['m'] if x[0] in used]+ghosts
    return g

def escher_geometry(m):
    """Original Escher drawing (pixel space): reactions as bezier edges between metabolite circles; None if the raw map is unavailable."""
    if m['source']!='escher':return None
    name=m['name'];path=INPUTS/'raw/escher'/name.split('.')[0]/(name+'.json')
    if not path.exists():return None
    if path not in _RAW:_RAW[path]=json.loads(path.read_text())[1]
    raw=_RAW[path];nodes=raw['nodes'];mets={};segs=[];labs=[]
    for r in m['reactions']:
        rr=raw['reactions'].get(str(r['key']).split(':',1)[-1])
        if not rr:return None
        coef={x['bigg_id']:x['coefficient'] for x in rr['metabolites']};rk=r['key']
        for sg in rr['segments'].values():
            a=nodes[sg['from_node_id']];b=nodes[sg['to_node_id']];ar=[]
            for n in(a,b):
                ar.append(int(n['node_type']=='metabolite' and (coef.get(n['bigg_id'],0)>0 or rr['reversibility'])))
                if n['node_type']=='metabolite':
                    nid=sg['from_node_id'] if n is a else sg['to_node_id'];mets[nid]=[nid,n['bigg_id'],round(n['x'],1),round(n['y'],1),int(bool(n.get('node_is_primary'))),round(n.get('label_x',n['x']),1),round(n.get('label_y',n['y']),1)]
            b1=sg.get('b1');b2=sg.get('b2')
            segs.append([round(a['x'],1),round(a['y'],1),round(b['x'],1),round(b['y'],1)]+([round(b1['x'],1),round(b1['y'],1),round(b2['x'],1),round(b2['y'],1)] if b1 and b2 else [None]*4)+[ar[0],ar[1],rk])
        labs.append([rk,round(rr['label_x'],1),round(rr['label_y'],1),rr['bigg_id']])
    return stub_far(dict(m=list(mets.values()),s=segs,l=labs),m)

SHARD=200
ZH={'Shape':'形状','chain':'链','cycle':'环','fork':'叉','Appendix reactions':'附属反应','Parent map':'父图','Sibling units':'同图兄弟单元','Geometry-only unit: the source map has no reaction identity, so the unit was cut by drawn shape only. Each drawn interaction between two metabolite glyphs is a pseudo-reaction (route waypoint, rows marked in the Nodes table); there are no reaction IDs or names.':'仅几何单元：来源图没有反应标识，所以该单元只按画出的形状切分。两个代谢物图元之间每条画出的交互视为一个伪反应（路由转折点，在“节点”表中标出）；没有反应ID或名称。','pseudo-reaction (drawn interaction)':'伪反应（画出的交互）','Module detail':'模块详情','Back to catalog':'返回图录','Loading…':'加载中…','No module selected.':'未选择模块。','Could not load module data.':'无法加载模块数据。',
'This module is not a variant representative; only representatives have detail pages (see display_ids.json).':'该模块不是变体代表，只有变体代表才有详情页（见 display_ids.json）。',
'Layout':'布局','Fit':'适应窗口','Zoom in':'放大','Zoom out':'缩小','Labels':'标签','Reaction':'反应','Metabolite':'代谢物','derived':'推导','observed':'观测',
'Drawn like the expert maps: circles = metabolites, lines with arrowheads = reactions (name beside the line). Pale = derived position. Drag to pan, wheel to zoom, hover to link graph and tables.':'画法与专家图一致：圆点 = 代谢物，带箭头的线 = 反应（反应名写在线旁）。浅色 = 推导位置。拖动平移，滚轮缩放，悬停可联动图与表格。',
'Context':'背景信息','Module code':'模块编号','Family':'家族','Variant':'变体（画法）','Source map':'来源图','Species':'物种','Licence':'许可','Quality':'质量','Size':'规模','Instances':'实例数',
'Same connectivity shape: {a} instances · {b} layout variants':'同连接形状 {a} 个实例 · {b} 种画法','Variant {a} of {b}':'画法 {a} / {b}','{a} reactions / route steps':'{a} 个反应／路径步','{a} instance(s) share this layout':'{a} 个实例共用此布局',
'Other layouts with the same connectivity shape':'同连接形状的其他画法','Previous variant':'上一个画法','Next variant':'下一个画法','Technical details (original IDs)':'技术细节（原始ID）',
'Reactions':'反应','Nodes':'节点','Ports and attachments':'端口与附属物','Code':'编号','Reaction ID':'反应ID','Name':'名称','Equation':'方程','Reversible':'可逆','Position':'位置','Compartments':'区室','References':'引用','yes':'是','no':'否',
'ID only':'仅ID','Role':'角色','Entity ID':'实体ID','Label / name':'标签／名称','Compartment':'区室','Port':'端口','Metabolite':'代谢物','Scope':'范围','Outside reactions':'外部反应数','Outside nodes':'外部节点数',
'Attachment glyph IDs (no names stored):':'附属物图元ID（未存储名称）：','none':'无','Names shown are exactly those stored in the source payload; where only an identifier exists it is shown and marked ID only.':'所示名称均为来源数据中保存的名称；只有标识符时直接显示标识符并标注“仅ID”。',
'role':'角色',"Redundancy": "冗余情况", "This layout is an independent structure in {a} of {b} instances.": "该布局在 {a} / {b} 个实例中是独立结构。", "Windows contained in this structure": "本结构包含的重叠窗口", "independent": "独立结构", "independent structure: not contained in a larger module of the same source map": "独立结构：不被同一来源图中更大的模块包含", "contained in a larger module of the same source map (overlapping window)": "被同一来源图中更大的模块包含（重叠窗口）", "contained in {c}": "包含于 {c}","Other drawings of the same reactions": "同一组反应的其他画法", "geometry only: no reaction identity, so the same reaction group cannot be determined": "仅几何，无反应身份，无法判定同一组反应", "No reaction group could be formed for this module (missing or duplicate reaction identities).": "该模块无法形成反应组（反应身份缺失或重复）。", "{m} modules in {k} source maps": "{m} 个模块，分布在 {k} 张来源图", "this reaction group spans {n} connectivity shapes (families)": "该反应组跨 {n} 种连接形状（家族）", "Reactions of this group": "本组反应", "No other map draws exactly these reactions": "没有其他地图画过同一组反应", "Drawn in other source maps": "在其他来源图中的画法", "Drawn in the same source map": "在同一张来源图中的画法",'Size tier':'规模分级','micro · 1 reaction':'微型 · 1 个反应','small · 2–3 reactions':'小型 · 2–3 个反应','structure · ≥4 reactions':'结构单元 · ≥4 个反应','micro: stored for provenance and statistics, not used for retrieval':'微型：仅为来源记录与统计保存，不参与检索','small: used only to fill gaps':'小型：只用于补缺','structure: reusable local layout':'结构单元：可复用的局部布局','+{n} more variants: show this family in the catalog':'另有 {n} 种画法：在图录中查看该家族','Only layout':'唯一画法','Same connectivity shape: {a} instances · only layout':'同连接形状 {a} 个实例 · 唯一画法','Line style':'线型','The module has {n} separate parts spread across the source map; Fit shows the largest part, Show all shows everything.':'该模块在来源图上分成 {n} 个彼此分离的部分；“适应窗口”显示最大的部分，“显示全部”显示所有部分。','Show all':'显示全部','Unified straight':'统一直线','Source style':'原始线型','Curved lines exist only in Escher maps and are not used by the layout algorithm; Unified straight draws every reaction path as straight segments.':'曲线只存在于 Escher 地图中，布局算法并不使用；“统一直线”把每条反应路径都画成直线段。','Connected to a metabolite far outside this module (stub only)':'连接到本模块之外很远的代谢物（仅画出短线头）','Route waypoint':'路径转折点','Failed to load data':'数据加载失败','Hint: keep the shards folder next to module.html, and do not move or rename files in results/catalog.':'提示：请保持 shards 文件夹与 module.html 在同一目录，不要移动或重命名 results/catalog 中的文件。','Timed out after 10 s while loading the module data.':'加载模块数据超过 10 秒仍未完成。','Reaction lines bend at a stored reaction point (placed at the reaction label position), which is not part of the original map drawing.':'反应线经过一个存储的反应点（放在反应标签位置），该点并不是原图中画出的内容。','Many nodes: labels are hidden by default; tick Labels to show them.':'节点较多：默认隐藏标签，勾选“标签”可显示。'}

DETAIL_JS='''<script>
(function(){
var $=function(i){return document.getElementById(i)},app=$('app');
function esc(x){return String(x).replace(/[&<>"]/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]})}
var id=(location.search.match(/[?&]id=([MU]\\d+)/)||[])[1];
if(!id){app.textContent=tr('No module selected.');return}
var k=Math.floor(parseInt(id.slice(1),10)/%(shard)d),R=null,D=[],fam=[];
function load(src,cb){var s=document.createElement('script');s.src=src;s.onload=cb;s.onerror=function(){if(src===shard){clearTimeout(timer);fail(tr('Could not load module data.'))}};document.head.appendChild(s)}
var shard='shards/s'+k+'.js',done=false;
function fail(msg){if(done&&R)return;app.innerHTML='<div class="note err"><b>'+tr('Failed to load data')+'</b><br>'+msg+'<br><code>'+shard+'</code><br>'+tr('Hint: keep the shards folder next to module.html, and do not move or rename files in results/catalog.')+'</div>'}
var timer=setTimeout(function(){if(!R)fail(tr('Timed out after 10 s while loading the module data.'))},10000);
load(shard,function(){done=true;clearTimeout(timer);go()});
function th(s){var p=s.split('|'),n=p[0],o='<svg class="th" viewBox="-3 -3 66 66" aria-hidden="true">',i,pts=[];
for(i=0;i<n.length;i+=3)pts.push([n.charCodeAt(i)-48,n.charCodeAt(i+1)-48,n.charCodeAt(i+2)-48]);
if(p[1])p[1].split(',').forEach(function(e){var a=e.split('-'),u=pts[a[0]],v=pts[a[1]];if(u&&v)o+='<line class="edge" x1="'+u[0]+'" y1="'+u[1]+'" x2="'+v[0]+'" y2="'+v[1]+'"/>'});
pts.forEach(function(q){if(q[2]&1)return;o+='<circle class="nm" cx="'+q[0]+'" cy="'+q[1]+'" r="'+(pts.length>60?2.2:3)+'"'+((q[2]&2)?' opacity=".4"':'')+'/>'});return o+'</svg>'}
function go(){try{go2()}catch(e){fail(esc(String(e&&e.message||e)))}}
function go2(){R=(window.__S&&__S[k]||{})[id];D=window.CATALOG||[];
 if(!R){app.textContent=tr('This module is not a variant representative; only representatives have detail pages (see display_ids.json).');return}
 var b=$('back');if(document.referrer&&/\\/catalog(_units)?\\//.test(document.referrer))b.onclick=function(e){e.preventDefault();history.back()};
 document.title=R.name+' · Open Meta Maps';document.querySelector('main h1').innerHTML='<span class="dn">'+esc(R.name)+'</span>';
 render();document.addEventListener('langchange',function(){try{render()}catch(e){fail(esc(String(e&&e.message||e)))}});
 if(window.CATALOG)famReady();else load('search-index.js',famReady)}
function famReady(){D=window.CATALOG||[];fam=D.filter(function(d){return d[9]===R.fid}).sort(function(a,b){return a[15]-b[15]});famDone=true;fillFam();fillGroup();fillRed()}
var ent={},rx={},nodeIdx={},GEO=null,famDone=false,tok=0,TB=[],LS='s';try{LS=localStorage.getItem('ls')==='o'?'o':'s'}catch(e){}
function label(n){return n[6]||n[5]}
function geom(){var G=R.g,u=1,m=[],s=[],l=[],rm={},pts=[];
 if(G){G.m.forEach(function(x){var n=nodeIdx['g:'+x[0]];m.push(x.concat(['g:'+x[0],n?n[4]:1,x[1]]));rm[x[2]+','+x[3]]=(x[4]?17:10)});s=G.s;l=G.l}
 else{u=.35;R.n.forEach(function(n){if(!n[1])m.push([n[0],n[6]||n[5],n[2]*100,n[3]*100,1,n[2]*100+8*u,n[3]*100-8*u,n[0],n[4],n[5]]);else{var r=rx[n[0]];if(r)l.push([n[0],n[2]*100+6*u,n[3]*100-6*u,r[1]])}});
  R.e.forEach(function(e){var a=nodeIdx[e[0]],b=nodeIdx[e[1]];if(!a||!b)return;var r=a[1]?a:b,mt=a[1]?b:a,prod=e[2]==='product';s.push([r[2]*100,r[3]*100,mt[2]*100,mt[3]*100,null,null,null,null,0,prod?1:0,r[0]]);});
  m.forEach(function(x){rm[x[2]+','+x[3]]=17*u})}
 var bb=[1e9,1e9,-1e9,-1e9];function P(x,y){bb=[Math.min(bb[0],x),Math.min(bb[1],y),Math.max(bb[2],x),Math.max(bb[3],y)]}
 m.forEach(function(x){P(x[2]-20*u,x[3]-20*u);P(x[2]+20*u+60*u,x[3]+20*u)});s.forEach(function(x){P(x[0],x[1]);P(x[2],x[3])});l.forEach(function(x){P(x[1],x[2]-14*u);P(x[1]+60*u,x[2])});
 var par={};function f(a){while(par[a]!==undefined&&par[a]!==a)a=par[a];return a}function un(a,b){a=f(a);b=f(b);if(a!==b)par[a]=b}
 s.forEach(function(x){var a=x[0]+','+x[1],b=x[2]+','+x[3];if(par[a]===undefined)par[a]=a;if(par[b]===undefined)par[b]=b;un(a,b)});
 var cnt={};s.forEach(function(x){var r=f(x[0]+','+x[1]);cnt[r]=(cnt[r]||0)+1});var best=null;Object.keys(cnt).forEach(function(r){if(best===null||cnt[r]>cnt[best])best=r});
 var mb=[1e9,1e9,-1e9,-1e9];s.forEach(function(x){if(f(x[0]+','+x[1])===best){mb=[Math.min(mb[0],x[0],x[2]),Math.min(mb[1],x[1],x[3]),Math.max(mb[2],x[0],x[2]),Math.max(mb[3],x[1],x[3])]}});
 var pad=40*u;mb=[mb[0]-pad,mb[1]-pad,mb[2]+pad+60*u,mb[3]+pad];var parts=Object.keys(cnt).length,main=(mb[2]-mb[0])*(mb[3]-mb[1])*4<(bb[2]-bb[0])*(bb[3]-bb[1])&&parts>1;
 return {m:m,s:s,l:l,u:u,rm:rm,bb:bb,mb:main?mb:bb,parts:parts,sparse:main}}
function buildSvg(GE,curved){var u=GE.u,svg='<svg id="gv" class="gv" preserveAspectRatio="xMidYMid meet">';
 GE.s.forEach(function(s){var cv=curved&&s[4]!==null,d=!cv?'M'+s[0]+','+s[1]+'L'+s[2]+','+s[3]:'M'+s[0]+','+s[1]+'C'+s[4]+','+s[5]+' '+s[6]+','+s[7]+' '+s[2]+','+s[3];
  svg+='<path class="redge'+(s[11]?' far':'')+'" data-k="'+esc(s[10])+'" d="'+d+'" style="stroke-width:'+3*u+'"/>';
  [[0,s[0],s[1],!cv?s[2]:s[4],!cv?s[3]:s[5]],[1,s[2],s[3],!cv?s[0]:s[6],!cv?s[1]:s[7]]].forEach(function(e){if(s[8+e[0]]){var dx=e[1]-e[3],dy=e[2]-e[4],l=Math.sqrt(dx*dx+dy*dy)||1;dx/=l;dy/=l;var rad=GE.rm[e[1]+','+e[2]]||0,tx=e[1]-dx*(rad+2*u),ty=e[2]-dy*(rad+2*u),nx=-dy,ny=dx;
   svg+='<polygon class="arrow" points="'+tx+','+ty+' '+(tx-dx*19*u+nx*7*u)+','+(ty-dy*19*u+ny*7*u)+' '+(tx-dx*19*u-nx*7*u)+','+(ty-dy*19*u-ny*7*u)+'"/>'}})});
 GE.l.forEach(function(l){svg+='<text class="lab rl" data-k="'+esc(l[0])+'" x="'+l[1]+'" y="'+l[2]+'" font-size="'+14*u+'">'+esc(l[3])+'</text>'});
 GE.m.forEach(function(m){var far=String(m[0]).slice(-1)==='~';svg+='<g data-k="'+esc(m[7])+'" opacity="'+(m[8]?1:.45)+'"><circle class="nm'+(far?' far':'')+'" cx="'+m[2]+'" cy="'+m[3]+'" r="'+(m[4]?17:10)*u+'"><title>'+esc(far?tr('Connected to a metabolite far outside this module (stub only)'):m[9])+'</title></circle><text class="lab" x="'+m[5]+'" y="'+m[6]+'" font-size="'+(m[4]?15:12)*u+'">'+esc(m[1])+(far?' ↗':'')+'</text></g>'});
 svg+='</svg>';
 return svg}
var TIER=[['micro','微型','1 reaction','1 个反应','micro: stored for provenance and statistics, not used for retrieval','微型：仅为来源记录与统计保存，不参与检索'],['small','小型','2–3 reactions','2–3 个反应','small: used only to fill gaps','小型：只用于补缺'],['structure','结构单元','≥4 reactions','≥4 个反应','structure: reusable local layout','结构单元：可复用的局部布局']];
function tierb(n){var t=TIER[n<=1?0:(n<=3?1:2)];return '<span class="tag tier" title="'+tr(t[4])+'">'+tr(t[0]+' · '+t[2])+'</span>'}
function redb(){var r=R.rd;if(!r)return '—';var h=r.mx>0?'<span class="tag ind" title="'+tr('independent structure: not contained in a larger module of the same source map')+'">'+tr('independent')+'</span>':'<a class="tag cont" href="module.html?id='+r.link+'" title="'+tr('contained in a larger module of the same source map (overlapping window)')+'">'+tr('contained in {c}',{c:r.cont})+'</a>';
 if(r.mx>0&&!r.rep_max)h+=' <span class="meta">'+tr('This layout is an independent structure in {a} of {b} instances.',{a:r.mx,b:r.n})+'</span>';return h}
function redPanel(){var r=R.rd;if(!r||!r.contained.length)return '';return '<section class="panel" id="redBox"><details><summary><b>'+tr('Windows contained in this structure')+' ('+r.contained.length+')</b></summary><div id="redList"><p class="meta">'+tr('Loading…')+'</p></div></details></section>'}
function fillRed(){var box=$('redList');if(!box||!famDone)return;var byM={};D.forEach(function(x){byM[x[21]]=x});box.innerHTML='<div class="strip">'+R.rd.contained.map(function(c){var x=byM[c[2]];return '<a class="sv" href="module.html?id='+c[2]+'">'+(x?th(x[17]):'')+'<span>'+c[0]+' · '+c[1]+'</span></a>'}).join('')+'</div>'}
function render(){
 GEO=null;TB=[];var N=R.n,i;ent={};rx={};nodeIdx={};
 N.forEach(function(n,i){nodeIdx[n[0]]=n;if(!ent[n[5]]&&n[6])ent[n[5]]=n[6]});R.r.forEach(function(r){rx[r[0]]=r});
 var GE=geom(),u=GE.u,svg=buildSvg(GE,LS==='o');
 var kv=function(a,b){return '<dt>'+tr(a)+'</dt><dd>'+b+'</dd>'};
 var ctx='<dl class="kv ctx">'+kv('Module code','<code>'+id+'</code>')+kv('Family',R.fc+' · '+(R.fvar>1?tr('Same connectivity shape: {a} instances · {b} layout variants',{a:R.finst,b:R.fvar}):tr('Same connectivity shape: {a} instances · only layout',{a:R.finst})))+kv('Variant',R.fc+'-v'+R.vno+' · '+(R.fvar>1?tr('Variant {a} of {b}',{a:R.vno,b:R.fvar}):tr('Only layout')))+kv('Source map','<span class="dn">'+esc(R.map)+'</span>')+kv('Species','<span class="dn">'+esc(R.species)+'</span>')+kv('Source',esc(R.source))+kv('Quality',esc(R.quality))+kv('Licence',esc(R.license))+(R.u?'':kv('Redundancy',redb())+kv('Size tier',tierb(R.size)))+(R.u?kv('Shape',tr(R.u.shape))+kv('Appendix reactions',String(R.u.appendix))+kv('Parent map','<span class="dn">'+esc(R.u.parent)+'</span>')+kv('Sibling units',String(R.u.siblings)):'')+kv('Size',tr('{a} reactions / route steps',{a:R.size})+' · '+tr('{a} instance(s) share this layout',{a:R.count}))+'</dl>'
  +'<details><summary>'+tr('Technical details (original IDs)')+'</summary><pre>'+esc(JSON.stringify({display_module:id,module_id:R.mid,display_family:R.fc,family_id:R.fid,display_variant:R.fc+'-v'+R.vno,variant_id:R.vid,source_map_id:R.mapid,source_sha256:R.sha},null,2))+'</pre></details>';
 var eq=function(r){var s=[],p=[];r[4].forEach(function(c){var nm=ent[c[0]]||c[0],v=Math.abs(c[1]),t=(v!==1?(Math.round(v*100)/100)+' ':'')+nm;(c[1]<0?s:p).push(t)});return '<span class="dn">'+esc(s.join(' + '))+' '+(r[3]?'⇌':'→')+' '+esc(p.join(' + '))+'</span>'};
 var idonly=' <span class="muted">('+tr('ID only')+')</span>';
 var rt=R.r.map(function(r,j){var nn=nodeIdx[r[0]],nm=r[2]&&r[2]!==r[1]?'<span class="dn">'+esc(r[2])+'</span>':'<span class="muted">—</span>'+idonly;
  var refs=Object.keys(r[5]||{}).map(function(a){return esc(a)+': '+esc([].concat(r[5][a]).join(', '))}).join('; ')||'—';
  return '<tr data-k="'+esc(r[0])+'"><td><code>'+id+'-R'+(j+1)+'</code></td><td>'+esc(r[1])+'</td><td>'+nm+'</td><td>'+eq(r)+'</td><td>'+tr(r[3]?'yes':'no')+'</td><td>'+(nn?tr(nn[4]?'observed':'derived'):'—')+'</td><td>'+esc((r[6]||[]).join(', ')||'—')+'</td><td>'+refs+'</td></tr>'});
 var nt=N.map(function(n,j){var e5=String(n[5]),cp=e5.indexOf('@')>=0?e5.split('@')[1]:'';
  return '<tr data-k="'+esc(n[0])+'"><td>'+(j+1)+'</td><td>'+tr(n[1]?'Reaction':'Metabolite')+'</td><td>'+(R.u&&n[1]&&String(n[5]).indexOf('route:')===0?'<span class="muted">'+tr('pseudo-reaction (drawn interaction)')+'</span>':esc(n[5]))+'</td><td>'+(n[6]?'<span class="dn">'+esc(n[6])+'</span>':'<span class="muted">—</span>'+idonly)+'</td><td>'+tr(n[4]?'observed':'derived')+'</td><td>'+esc(cp||'—')+'</td></tr>'});
 var pt=R.p.map(function(p){return '<tr data-k="'+esc(p[0]||'')+'"><td>'+esc(p[1]||p[0]||'—')+(p[1]?'':idonly)+'</td><td>'+esc(p[2]||'—')+'</td><td>'+p[3]+'</td></tr>'});
 var tb=function(h,rows){var i=TB.length;TB.push(rows);return '<div class="tablewrap"><table><thead><tr>'+h.map(function(x){return '<th>'+tr(x)+'</th>'}).join('')+'</tr></thead><tbody id="tb'+i+'"></tbody></table></div>'};
 var big=(GE.m.length+GE.l.length)>300,checked=($('lb')?$('lb').checked:!big);
 app.innerHTML='<p class="note">'+tr('Drawn like the expert maps: circles = metabolites, lines with arrowheads = reactions (name beside the line). Pale = derived position. Drag to pan, wheel to zoom, hover to link graph and tables.')+(R.g?'':' '+tr('Reaction lines bend at a stored reaction point (placed at the reaction label position), which is not part of the original map drawing.'))+(GE.sparse?' '+tr('The module has {n} separate parts spread across the source map; Fit shows the largest part, Show all shows everything.',{n:GE.parts}):'')+(big?' '+tr('Many nodes: labels are hidden by default; tick Labels to show them.'):'')+'</p>'
  +'<div class="viewer"><div class="vbar"><button type="button" id="vf">'+tr('Fit')+'</button><button type="button" id="va">'+tr('Show all')+'</button><button type="button" id="vi" aria-label="'+tr('Zoom in')+'">+</button><button type="button" id="vo" aria-label="'+tr('Zoom out')+'">−</button><label>'+tr('Line style')+' <select id="ls"><option value="s"'+(LS==='s'?' selected':'')+'>'+tr('Unified straight')+'</option><option value="o"'+(LS==='o'?' selected':'')+'>'+tr('Source style')+'</option></select></label><label><input type="checkbox" id="lb"'+(checked?' checked':'')+'> '+tr('Labels')+'</label></div><p class="meta lsn">'+tr('Curved lines exist only in Escher maps and are not used by the layout algorithm; Unified straight draws every reaction path as straight segments.')+'</p><div class="vbox'+(checked?'':' nolab')+'" id="vbox">'+svg+'</div></div>'
  +'<div class="detail-grid'+(R.fvar>1?'':' one')+'"><section class="panel"><h2>'+tr('Context')+'</h2>'+ctx+'</section>'+(R.fvar>1?'<section class="panel"><h2>'+tr('Other layouts with the same connectivity shape')+'</h2><div id="famBox"><p class="meta">'+tr('Loading…')+'</p></div></section>':'')+'</div>'
  +grpPanel()+redPanel()
  +'<h2>'+tr('Reactions')+'</h2>'+(R.u&&!R.r.length?'<p class="note">'+tr('Geometry-only unit: the source map has no reaction identity, so the unit was cut by drawn shape only. Each drawn interaction between two metabolite glyphs is a pseudo-reaction (route waypoint, rows marked in the Nodes table); there are no reaction IDs or names.')+'</p>':'')+tb(['Code','Reaction ID','Name','Equation','Reversible','Position','Compartments','References'],rt)
  +'<h2>'+tr('Nodes')+'</h2>'+tb(['#','Role','Entity ID','Label / name','Position','Compartment'],nt)
  +'<h2>'+tr('Ports and attachments')+'</h2>'+(pt?tb(['Port','Role','Outside reactions'],pt):'<p class="meta">'+tr('none')+'</p>')
  +'<p class="meta">'+tr('Attachment glyph IDs (no names stored):')+' <code>'+(R.a.length?esc(R.a.join(', ')):tr('none'))+'</code></p>'
  +'<p class="meta">'+tr('Names shown are exactly those stored in the source payload; where only an identifier exists it is shown and marked ID only.')+'</p>';
 view();fillTables();fillFam();fillGroup();fillRed()}
function grpPanel(){if(R.u)return '';var gp=R.gp,h='<section class="panel" id="grpBox"><h2>'+tr('Other drawings of the same reactions')+'</h2>';
 if(R.quality==='geometry_only'||!R.r.length)return h+'<p class="meta">'+tr('geometry only: no reaction identity, so the same reaction group cannot be determined')+'</p></section>';
 if(!gp)return h+'<p class="meta">'+tr('No reaction group could be formed for this module (missing or duplicate reaction identities).')+'</p></section>';
 var ids=gp.reactions,short=ids.slice(0,6).join(', ')+(ids.length>6?' …':'');
 h+='<p class="meta"><code>'+esc(gp.code)+'</code> · '+tr('{m} modules in {k} source maps',{m:gp.nm,k:gp.nmaps})+(gp.fams>1?' · '+tr('this reaction group spans {n} connectivity shapes (families)',{n:gp.fams}):'')+'</p>'
  +'<details'+(ids.length<=6?' open':'')+'><summary>'+tr('Reactions of this group')+': <span class="dn">'+esc(short)+'</span></summary><pre>'+esc(ids.join('\\n'))+'</pre></details>';
 if(!gp.p.length)return h+'<p class="meta">'+tr('No other map draws exactly these reactions')+'</p></section>';
 return h+'<div id="grpList"><p class="meta">'+tr('Loading…')+'</p></div></section>'}
function fillGroup(){var box=$('grpList');if(!box||!famDone)return;var gp=R.gp,byM={};D.forEach(function(x){byM[x[21]]=x});
 var other=gp.p.filter(function(p){return !p[4]}),same=gp.p.filter(function(p){return p[4]});
 function block(list,title){if(!list.length)return '';var g={};list.forEach(function(p){var k=p[1]+'|'+p[2];(g[k]=g[k]||[]).push(p)});
  return '<h3>'+title+'</h3>'+Object.keys(g).map(function(k){var l=g[k],p0=l[0];return '<div class="pgrp"><div class="pt"><span class="dn">'+esc(p0[1])+'</span> · <span class="dn">'+esc(p0[2])+'</span> · '+esc(p0[3])+'</div><div class="strip">'+l.map(function(p){var r=byM[p[0]];return '<a class="sv" href="module.html?id='+p[0]+'"'+(p[0]===id?' aria-current="true"':'')+'>'+(r?th(r[17]):'')+'<span>'+p[0]+'</span></a>'}).join('')+'</div></div>'}).join('')}
 box.innerHTML=block(other,tr('Drawn in other source maps'))+block(same,tr('Drawn in the same source map'))}
function fillTables(){var my=++tok;TB.forEach(function(rows,i){var el=$('tb'+i),p=0;(function step(){if(my!==tok||!el)return;el.insertAdjacentHTML('beforeend',rows.slice(p,p+60).join(''));p+=60;if(p<rows.length)setTimeout(step,0)})()})}
function fillFam(){var box=$('famBox');if(!box||!famDone)return;var cur=0;fam.forEach(function(d,j){if(d[21]===id)cur=j});var shown=fam.slice(0,12);if(cur>=12)shown.push(fam[cur]);var more=fam.length-shown.length,strip=shown.map(function(d){return '<a class="sv" href="module.html?id='+d[21]+'"'+(d[21]===id?' aria-current="true"':'')+'>'+th(d[17])+'<span>'+d[14]+'-v'+d[15]+'</span></a>'}).join(''),prev=fam[cur-1],next=fam[cur+1];
 box.innerHTML='<div class="pager2">'+(prev?'<a href="module.html?id='+prev[21]+'">← '+tr('Previous variant')+'</a>':'<span></span>')+(next?'<a href="module.html?id='+next[21]+'">'+tr('Next variant')+' →</a>':'<span></span>')+'</div><div class="strip">'+strip+'</div>'+(more>0?'<p class="meta"><a href="index.html#f='+R.fid+'">'+tr('+{n} more variants: show this family in the catalog',{n:more})+'</a></p>':'')}
function view(keep){var box=$('vbox'),sv=$('gv'),N=R.n,vb,fit;
 function setvb(){sv.setAttribute('viewBox',vb.join(' '))}
 function doFit(all){var G0=geom(),bb=all===true?G0.bb:G0.mb,x0=bb[0],y0=bb[1],x1=bb[2],y1=bb[3];
  var ext=Math.max(x1-x0,y1-y0,60),pad=ext*.04+10,w=x1-x0+2*pad,h=y1-y0+2*pad,a=(box.clientHeight||420)/(box.clientWidth||800);if(h/w<a)h=w*a;else w=h/a;vb=[(x0+x1)/2-w/2,(y0+y1)/2-h/2,w,h];setvb()}
 function zoom(f,cx,cy){var px=cx===undefined?vb[0]+vb[2]/2:cx,py=cy===undefined?vb[1]+vb[3]/2:cy;vb=[px-(px-vb[0])*f,py-(py-vb[1])*f,vb[2]*f,vb[3]*f];setvb()}
 function pt(e){var r=sv.getBoundingClientRect();return [vb[0]+(e.clientX-r.left)/r.width*vb[2],vb[1]+(e.clientY-r.top)/r.height*vb[3]]}
 if(keep){vb=keep.split(' ').map(Number);setvb()}else doFit();$('ls').onchange=function(){LS=this.value;try{localStorage.setItem('ls',LS)}catch(e){}var kv=sv.getAttribute('viewBox');box.innerHTML=buildSvg(geom(),LS==='o');view(kv)};$('vf').onclick=function(){doFit()};$('va').onclick=function(){doFit(true)};$('vi').onclick=function(){zoom(.7)};$('vo').onclick=function(){zoom(1.4)};
 $('lb').onchange=function(){box.classList.toggle('nolab',!this.checked)};
 sv.addEventListener('wheel',function(e){e.preventDefault();var p=pt(e);zoom(e.deltaY<0?.85:1.18,p[0],p[1])},{passive:false});
 var drag=null;sv.addEventListener('pointerdown',function(e){drag={x:e.clientX,y:e.clientY,v:vb.slice()};sv.setPointerCapture(e.pointerId)});
 sv.addEventListener('pointermove',function(e){if(!drag)return;var r=sv.getBoundingClientRect();vb=[drag.v[0]-(e.clientX-drag.x)/r.width*drag.v[2],drag.v[1]-(e.clientY-drag.y)/r.height*drag.v[3],drag.v[2],drag.v[3]];setvb()});
 sv.addEventListener('pointerup',function(){drag=null});
 window.onresize=function(){doFit()}}
document.addEventListener('mouseover',function(e){var t=e.target.closest&&e.target.closest('[data-k]');document.querySelectorAll('.hl').forEach(function(x){x.classList.remove('hl')});if(!t)return;var key=t.getAttribute('data-k');document.querySelectorAll('[data-k]').forEach(function(x){if(x.getAttribute('data-k')===key)x.classList.add('hl')})});
})();
</script>''' % {'shard':SHARD}

def shard_key(code):return int(code[1:])//SHARD

def text(x):
    """Display string for stored identifiers: xref dicts become 'Database:ID', missing values ''."""
    if x is None:return ''
    if isinstance(x,dict):return ':'.join(str(x[k]) for k in ('Database','ID') if k in x) or json.dumps(x,sort_keys=True)
    return str(x)

def record(m,mc,fc,vno,finst,fvar,count,gp=None,rd=None,lib='windows'):
    extra=dict(u=dict(shape=m['shape_type'],appendix=len(m.get('appendix_reactions',[])),parent=m.get('parent_map',{}).get('name',m['name']),siblings=len(m.get('sibling_units',[])))) if lib=='units' else {}
    return dict(extra,gp=gp,rd=rd,name=m['name'],mid=m['id'],fid=m['core_family'],vid=m['variant'],mapid=m['map_id'],sha=m['source_sha256'],map=m['name'],species=m['species'],source=m['source'],quality=m['quality'],license=m['license'],size=m['size'],fc=fc,vno=vno,finst=finst,fvar=fvar,count=count,
     n=[[n['key'],int(n['role']=='reaction'),round(n['xy'][0],3),round(n['xy'][1],3),int(n['observed']),text(n['entity']),text(n.get('label'))] for n in m['nodes']],
     e=[[e['a'],e['b'],e.get('flow','')] for e in m['edges']],g=escher_geometry(m),
     r=[[x['key'],text(x['id']),text(x.get('name')),int(bool(x.get('reversible'))),[[text(k),v] for k,v in x['coefficients'].items()],x.get('references') or {},x.get('compartments') or []] for x in m['reactions']],
     p=[[x.get('key'),x.get('metabolite') or '',x.get('role',''),len(x.get('outside_reactions',x.get('outside_nodes',[])))] for x in m['ports']],a=m['attachments'])

def write(dest,items,lib='windows'):
    """items: iterable of (module, module_code, family_code, variant_no, family_instances, family_variants, instances, reaction_group_or_None)."""
    sd=dest/'shards';sd.mkdir(exist_ok=True);groups={}
    for m,mc,fc,vno,fi,fv,c,gp,rd in items:groups.setdefault(shard_key(mc),{})[mc]=record(m,mc,fc,vno,fi,fv,c,gp,rd,lib)
    for p in sd.glob('s*.js'):p.unlink()
    for k,v in groups.items():(sd/f's{k}.js').write_text('(window.__S=window.__S||{})['+str(k)+']='+json.dumps(v,ensure_ascii=False,separators=(',',':'))+';')
    body=('<div class="toolbar"><a id="back" href="index.html">← <span>Back to catalog</span></a></div><div id="app"><p>Loading…</p></div>')
    ui.write(dest/'module.html',('Module detail','模块详情'),body,section='catalog',zh=ZH,script=DETAIL_JS,crumbs=[('Catalog','index.html'),(('Module detail','模块详情'),None)])
    return len(groups)
