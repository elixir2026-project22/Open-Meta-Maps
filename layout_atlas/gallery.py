"""Paginated actual-instance catalog and portable compressed JSON exchange."""
import gzip,html
from core import *
import database,ui,gallery_detail,gallery_method

PER_PAGE=80
ZH={"Independent structures only": "只显示独立结构", "Independent structure = a module not contained in a larger module of the same source map; the others are overlapping windows": "独立结构 = 不被同一来源图中更大模块包含的模块；其余是重叠窗口", "independent": "独立结构", "independent structure: not contained in a larger module of the same source map": "独立结构：不被同一来源图中更大的模块包含", "contained in a larger module of the same source map (overlapping window)": "被同一来源图中更大的模块包含（重叠窗口）", "contained in {c}": "包含于 {c}", "{n} contained in a larger module": "{n} 个被更大模块包含", "{n} by size tier": "{n} 个因规模分级", "hidden": "已隐藏", "counts follow the filters above; the summary cards count all modules": "计数随上方筛选变化；汇总卡片统计全部模块", "In the current selection: {f} families, {v} layout variants (over all modules: {F} families, {V} layout variants)": "当前选择中：{f} 个家族、{v} 个布局变体（全部模块：{F} 个家族、{V} 个布局变体）","reaction set drawn in ≥2 maps": "同一组反应在 ≥2 张图中出现", "Modules whose exact set of reactions is also drawn in at least one other source map (lower bound)": "完全相同的反应集合在至少另一张来源图中也被画过（下限值）", "same reactions drawn in more than one source map": "同一组反应在多张来源图中被画过", "same reactions · {n} maps": "同组反应 · {n} 张图","Size tier": "规模分级", "Show all": "显示全部", "Showing {a} of {b} module instances": "显示 {a} / {b} 个模块实例", " (micro hidden: {c})": "（已隐藏微型：{c}）", "counts follow the tier selection above; the summary cards count all tiers": "计数随上方分级选择变化；上方汇总卡片统计全部分级", "micro · 1 reaction": "微型 · 1 个反应", "small · 2–3 reactions": "小型 · 2–3 个反应", "structure · ≥4 reactions": "结构单元 · ≥4 个反应", "micro: stored for provenance and statistics, not used for retrieval": "微型：仅为来源记录与统计保存，不参与检索", "small: used only to fill gaps": "小型：只用于补缺", "structure: reusable local layout": "结构单元：可复用的局部布局","Only layout":"唯一画法","Same connectivity shape: {a} instances · only layout":"同连接形状 {a} 个实例 · 唯一画法","Open detail view":"打开详情页","module id":"模块ID","source map id":"来源图ID","Technical details (original IDs)":"技术细节（原始ID）","Module code":"模块编号","Same colour = same connectivity shape (family; reactions and metabolites can differ); variants = different layouts of that shape.": "相同颜色 = 相同连接形状（家族，反应和代谢物可以不同）；变体 = 该形状的不同画法。", "What do these terms mean?": "这些术语是什么意思？", "Show all families": "显示全部家族", "Family · variant": "家族 · 变体", "Family = modules with the same connectivity shape (reactions and metabolites can differ; only the connection pattern is the same). Variant = a different drawn layout of that shape.": "家族 = 连接形状相同的模块（反应和代谢物可以不同，只是连接模式相同）；变体（画法）= 该形状的不同布局画法。", "geometry only, no reaction IDs": "仅几何，无反应ID", "Click to show only this family": "点击仅显示该家族", "Variant {a} of {b}": "画法 {a} / {b}", "Same connectivity shape: {a} instances · {b} layout variants": "同连接形状 {a} 个实例 · {b} 种画法", "Showing one family only": "仅显示一个家族", "{a} of {b} layout variants match": "{a} / {b} 个布局变体匹配", " (first {n} shown)": "（显示前 {n} 个）","Layout": "布局缩略图", "Module and reactions": "模块与反应", "Family · variant": "家族 · 变体", "Order": "排序", "Group by family": "按家族分组", "Size (large first)": "规模（大者在前）", "Name A–Z": "名称 A–Z", "name, species, family or reaction id": "名称、物种、家族或反应ID", "name, species, family id or reaction id": "名称、物种、家族或反应ID","Shape":"形状","chain":"链","cycle":"环","fork":"叉","Appendix":"附属反应","shape: chain, cycle or fork":"形状：链、环或叉","Structure-unit reuse (experimental)":"结构单元复用（实验性）","Units never overlap, so there is no redundancy index; each unit belongs to exactly one source map.":"单元互不重叠，所以没有冗余索引；每个单元只属于一张来源图。","Expert metabolic-map structure units (non-overlapping, 3–12 reactions) cut once from Escher community maps, WikiPathways and Reactome. Experimental: the tool uses pathway-fragment reuse by default.":"从 Escher 社区地图、WikiPathways 和 Reactome 各切一次得到的专家代谢图结构单元（互不重叠，3–12 个反应）。实验性：工具默认使用通路片段复用。","Structure-unit catalog":"结构单元图录","Units by shape":"按形状统计的单元","{n} appendix reactions":"{n} 个附属反应","pseudo-reactions (drawn interactions), no reaction IDs":"伪反应（画出的交互），无反应ID"}
SRC_COLOR={'community':'#2d8f6f','wikipathways':'#3b77b8','reactome':'#a8559a'}

def rx_ids(m):return [r['id'] for r in m['reactions']]

def rx_summary(m,k=5):
 ids=rx_ids(m);s=', '.join(ids[:k]);return s+(f' +{len(ids)-k} more' if len(ids)>k else '')

def rx_html(m,k=5):
 ids=rx_ids(m)
 if not ids:return ui.bi('geometry only, no reaction IDs','仅几何，无反应ID')
 return '<span class="dn">'+html.escape(', '.join(ids[:k]))+'</span>'+(' '+ui.bi(f'+{len(ids)-k} more',f'另有 {len(ids)-k} 个') if len(ids)>k else '')

def map_html(m):
 n=m['name'];e=html.escape
 if m['source']=='community' and '.' in n:return ui.bi('Escher map','Escher 地图')+' '+e(n.split('.',1)[0])+' · <span class="dn">'+e(n.split('.',1)[1])+'</span>'
 if m['source']=='reactome':return 'Reactome '+e(m['lineage'].split(':',1)[-1])+' · <span class="dn">'+e(n)+'</span>'
 return e(m['source'] if m['source']!='wikipathways' else 'WikiPathways')+' · <span class="dn">'+e(n)+'</span>'

def map_line(m):
 n=m['name'];lin=m['lineage']
 if m['source']=='community' and '.' in n:return 'Escher map '+n.split('.',1)[0]+' · '+n.split('.',1)[1]
 if m['source']=='reactome':return 'Reactome '+lin.split(':',1)[-1]+' · '+n
 if m['source']=='wikipathways':return 'WikiPathways · '+n
 return m['source']+' · '+n

TIERS=[('micro','微型',1,'1 reaction','1 个反应','micro: stored for provenance and statistics, not used for retrieval','微型：仅为来源记录与统计保存，不参与检索'),('small','小型',2,'2–3 reactions','2–3 个反应','small: used only to fill gaps','小型：只用于补缺'),('structure','结构单元',4,'≥4 reactions','≥4 个反应','structure: reusable local layout','结构单元：可复用的局部布局')]

def tier(size):return 0 if size<=1 else (1 if size<=3 else 2)

def tier_badge(size):
 en,zh,_,ren,rzh,ten,tzh=TIERS[tier(size)]
 return f'<span class="tag tier t{tier(size)}" title="{html.escape(ten)}" data-title-zh="{html.escape(tzh)}">'+ui.bi(f'{en} · {ren}',f'{zh} · {rzh}')+'</span>'

CODES={};GROUPS={};RED={};CUR=['windows']
SHAPE_ZH={'chain':'链','cycle':'环','fork':'叉'};SHAPE_COLOR={'chain':'#3b77b8','cycle':'#2d8f6f','fork':'#b66b28'}

def is_units():return CUR[0]=='units'

def shape_tags(m):
 s=m.get('shape_type','');a=len(m.get('appendix_reactions',[]))
 h=f'<span class="tag" title="shape: chain, cycle or fork" data-title-zh="形状：链、环或叉">'+ui.bi(s,SHAPE_ZH.get(s,s))+'</span>'
 if a:h+='<span class="tag">'+ui.bi(f'+{a} appendix',f'+{a} 个附属反应')+'</span>'
 return h

def load_redundancy(mods,codes):
 """Read the derived module-redundancy index (read-only): per variant the independent-instance count and the container of its representative."""
 path=OUT/'module_redundancy.json.gz'
 if not path.exists():return {},{}
 with gzip.open(path,'rt') as f:red=json.load(f)
 by={m['id']:m for m in mods};rep_of={};n=defaultdict(int);mx=defaultdict(int)
 for m in mods:
  rep_of.setdefault(m['variant'],m['id']);n[m['variant']]+=1;mx[m['variant']]+=int(red[m['id']]['maximal'])
 info={}
 for v,rid in rep_of.items():
  r=red[rid];c=r['container'];cc=codes['m'][c] if c else '';link=codes['m'][rep_of[by[c]['variant']]] if c else ''
  info[rid]=dict(mx=mx[v],n=n[v],rep_max=int(r['maximal']),cont=cc,link=link)
 contained=defaultdict(list)
 for m in mods:
  c=red[m['id']]['container']
  if c and c in info:contained[c].append([codes['m'][m['id']],m['size'],codes['m'][rep_of[m['variant']]]])
 for rid in info:info[rid]['contained']=sorted(contained.get(rid,[]),key=lambda x:(-x[1],x[0]))
 return info,red

def red_badge(r):
 if not r:return ''
 if r['mx']>0:return '<span class="tag ind" title="independent structure: not contained in a larger module of the same source map" data-title-zh="独立结构：不被同一来源图中更大的模块包含">'+ui.bi('independent','独立结构')+'</span>'
 return '<a class="tag cont" href="module.html?id='+r['link']+'" title="contained in a larger module of the same source map (overlapping window)" data-title-zh="被同一来源图中更大的模块包含（重叠窗口）">'+ui.bi('contained in '+r['cont'],'包含于 '+r['cont'])+'</a>'

def group_badge(gp):
 if not gp or not gp['cross']:return ''
 n=gp['nmaps'];return '<span class="tag grp" title="same reactions drawn in more than one source map" data-title-zh="同一组反应在多张来源图中被画过">'+ui.bi(f'same reactions · {n} maps',f'同组反应 · {n} 张图')+'</span>'

def assign_codes(mods,families,fp='F',mp='M'):
 """Display-only codes: F0001.. by (instances desc, id); variants by (instances desc, id); M00001.. by (source, family, variant, id)."""
 fcode={f:fp+'%04d'%(i+1) for i,f in enumerate(sorted(families,key=lambda f:(-len(families[f]),f)))};vno={}
 for f,ins in families.items():
  c=Counter(m['variant'] for m in ins)
  for i,v in enumerate(sorted(c,key=lambda v:(-c[v],v)),1):vno[v]=i
 mcode={m['id']:mp+'%05d'%(i+1) for i,m in enumerate(sorted(mods,key=lambda m:(m['source'],fcode[m['core_family']],vno[m['variant']],m['id'])))}
 CODES.update(f=fcode,v=vno,m=mcode);return CODES

def load_groups(mods,codes):
 """Read the derived reaction-group index (read-only): per module its group code, canonical reactions and partner modules."""
 path=OUT/'reaction_groups.json.gz'
 if not path.exists():return {},{}
 with gzip.open(path,'rt') as f:rg=json.load(f)
 by={m['id']:m for m in mods};rep_of={}
 for m in mods:rep_of.setdefault(m['variant'],m['id'])
 order=sorted(rg['groups'],key=lambda x:(-x['source_maps'],-len(x['modules']),x['id']));gcode={x['id']:'G%04d'%(i+1) for i,x in enumerate(order)}
 out={};
 for x in rg['groups']:
  for mid in x['modules']:
   partners=[];seen=set()
   for o in sorted(x['modules']):
    if o==mid:continue
    om=by[o];rm=codes['m'][rep_of[om['variant']]];key=(rm,om['map_id'])
    if key in seen:continue
    seen.add(key);partners.append([rm,om['name'],om['species'],om['source'],int(om['map_id']==by[mid]['map_id'])])
   out[mid]=dict(code=gcode[x['id']],id=x['id'],reactions=x['reactions'],nm=len(x['modules']),nmaps=x['source_maps'],cross=int(x['cross_map']),fams=x['families'],p=partners)
 return out,gcode

def hue(family):return int(family[:4],16)%360

def thumb(m):
 ns=m['nodes'];xy=np.array([n['xy'] for n in ns],float);lo=xy.min(0);span=max((xy.max(0)-lo).max(),1e-9);q=np.rint((xy-lo)/span*60).astype(int);idx={n['key']:i for i,n in enumerate(ns)}
 nodes=''.join(chr(48+int(x))+chr(48+int(y))+chr(48+(n['role']=='reaction')+2*(not n['observed'])) for (x,y),n in zip(q,ns))
 edges=','.join(f"{idx[e['a']]}-{idx[e['b']]}" for e in m['edges'] if e['a'] in idx and e['b'] in idx)
 return nodes+'|'+edges

def bar(counts,colors=SRC_COLOR):
 total=sum(counts.values())
 segs=''.join(f'<i style="width:{100*v/total:.2f}%;background:{colors.get(k,"#888")}"></i>' for k,v in counts.items())
 legend=''.join(f'<em style="--c:{colors.get(k,"#888")}">{html.escape(k)} {v:,}</em>' for k,v in counts.items())
 return f'<div class="bar" role="img" aria-label="by source">{segs}</div><div class="legend">{legend}</div>'

def card(m,count,source_support,row,vno,finst=0,fvar=0):
 ns=m['nodes'];xy={n['key']:np.array(n['xy'])*55 for n in ns};a=np.array(list(xy.values()));lo=a.min(0)-15;size=a.max(0)-lo+15
 svg=[f'<svg viewBox="{lo[0]} {lo[1]} {size[0]} {size[1]}" role="img" aria-label="{html.escape(m["name"])}">']
 for e in m['edges']:
  p,q=xy[e['a']],xy[e['b']];svg.append(f'<path class="edge" d="M{p[0]},{p[1]}L{q[0]},{q[1]}"/>')
 for n in ns:
  p=xy[n['key']]
  if n['role']=='reaction':continue
  cls='nm';svg.append(f'<circle class="{cls}" cx="{p[0]}" cy="{p[1]}" r="3" opacity="{1 if n["observed"] else .4}"><title>{html.escape(str(n["entity"]))}: {n["position_origin"]}</title></circle>')
 svg.append('</svg>');e=html.escape
 fc=CODES['f'][m['core_family']];mc=CODES['m'][m['id']];vc=f'{fc}-v{vno}'
 text=' '.join([m['name'],m['species'],m['source'],m['core_family'],m['id'],fc,vc,mc]+rx_ids(m)).lower()
 return (f'<article class="card" id="m-{m["id"]}" data-src="{e(m["source"])}" data-sp="{e(m["species"])}" data-q="{e(m["quality"])}" data-t="{e(text)}">'
  f'<h2 class="dn"><a href="module.html?id={mc}">{e(m["name"])}</a></h2><p class="meta" title="source map id {e(m["map_id"])}" data-title-zh="来源图ID {e(m["map_id"])}">{map_html(m)}</p><a class="svgl" href="module.html?id={mc}" title="Open detail view" data-title-zh="打开详情页">{"".join(svg)}</a>'
  f'<p class="rx"><span class="fam" style="--h:{hue(m["core_family"])}" title="family id {e(m["core_family"])} · variant id {e(m["variant"])}" data-title-zh="家族ID {e(m["core_family"])} · 变体ID {e(m["variant"])}">{vc}</span> {rx_html(m)}</p>'
  f'<div class="tags"><span class="tag src-{e(m["source"])}">{e(m["source"])}</span><span class="tag">{e(m["species"])}</span><span class="tag">{e(m["quality"])}</span><span class="tag lic">{e(m["license"])}</span>{shape_tags(m) if is_units() else tier_badge(m["size"])+group_badge(GROUPS.get(m["id"]))+red_badge(RED.get(m["id"]))}</div>'
  f'<p class="meta">'+ui.bi(f'{m["size"]} reactions / route steps · {count} instance(s) · {source_support} source layout(s), not independent authors · {len(m["ports"])} ports · {len(m["attachments"])} attachments',f'{m["size"]} 个反应／路径步 · {count} 个实例 · {source_support} 个来源布局（并非独立作者） · {len(m["ports"])} 个端口 · {len(m["attachments"])} 个附件')+'</p>'
  f'<dl class="kv"><dt>Module code</dt><dd><code>{mc}</code></dd><dt>Family</dt><dd>{fc} · '+(ui.bi(f'Same connectivity shape: {finst} instances · {fvar} layout variants',f'同连接形状 {finst} 个实例 · {fvar} 种画法') if fvar>1 else ui.bi(f'Same connectivity shape: {finst} instances · only layout',f'同连接形状 {finst} 个实例 · 唯一画法'))+f'</dd><dt>Variant</dt><dd>'+(ui.bi(f'Variant {vno} of {fvar}',f'画法 {vno} / {fvar}') if fvar>1 else ui.bi('Only layout','唯一画法'))+f' ({vc})</dd><dt>Source map</dt><dd>{map_html(m)}</dd><dt>Licence</dt><dd>{e(m["license"])}</dd></dl>'
  f'<details><summary>Technical details (original IDs)</summary><pre>{e(json.dumps(dict(row,display_module=mc,display_family=fc,display_variant=vc),ensure_ascii=False,indent=2))}</pre>'
  f'<p>'+ui.bi('Full coordinates, paths and roles:','完整坐标、路径与角色：')+f' <code>modules.json.gz</code> '+ui.bi(f'record {m["id"]}.',f'中的记录 {m["id"]}。')+' '+ui.bi('Automatic evidence-based validation; no human certification claimed.','基于证据的自动验证；不声称经过人工认证。')+'</p></details></article>')

PAGE_JS='''<script>
(function(){var q=document.getElementById('q'),s=document.getElementById('src'),k=document.getElementById('cnt'),cs=[].slice.call(document.querySelectorAll('.card'));
function f(){var t=q.value.trim().toLowerCase(),v=s.value,n=0;cs.forEach(function(c){var ok=(!v||c.dataset.src===v)&&(!t||c.dataset.t.indexOf(t)>=0);c.style.display=ok?'':'none';if(ok)n++});k.textContent=tr('{a} of {b} cards on this page shown',{a:n,b:cs.length})}document.addEventListener('langchange',f);
q.oninput=s.onchange=f;f();
if(location.hash){var h=document.querySelector(location.hash);if(h){h.style.outline='2px solid var(--accent)'}}})();
</script>'''

INDEX_JS='''<script src="search-index.js"></script><script>
(function(){var D=window.CATALOG||[],$=function(i){return document.getElementById(i)},out=$('hits'),cnt=$('cnt'),LIM=120,F='';
function esc(x){return String(x).replace(/[&<>"]/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]})}
function maptxt(s){var m=s.match(/^Escher map (\\S+) · ([\\s\\S]*)$/);return m?tr('Escher map')+' '+esc(m[1])+' · <span class="dn">'+esc(m[2])+'</span>':'<span class="dn">'+esc(s)+'</span>'}
function th(s){var p=s.split('|'),n=p[0],o='<svg class="th" viewBox="-3 -3 66 66" aria-hidden="true">',i,pts=[];
for(i=0;i<n.length;i+=3)pts.push([n.charCodeAt(i)-48,n.charCodeAt(i+1)-48,n.charCodeAt(i+2)-48]);
if(p[1])p[1].split(',').forEach(function(e){var a=e.split('-'),u=pts[a[0]],v=pts[a[1]];if(u&&v)o+='<line class="edge" x1="'+u[0]+'" y1="'+u[1]+'" x2="'+v[0]+'" y2="'+v[1]+'"/>'});
pts.forEach(function(q){if(q[2]&1)return;o+='<circle class="nm" cx="'+q[0]+'" cy="'+q[1]+'" r="'+(pts.length>60?2.2:3)+'"'+((q[2]&2)?' opacity=".4"':'')+'/>'});return o+'</svg>'}
var TT=[['micro · 1 reaction','micro: stored for provenance and statistics, not used for retrieval'],['small · 2–3 reactions','small: used only to fill gaps'],['structure · ≥4 reactions','structure: reusable local layout']];var TOTAL=%(total)d,MICRO=%(micro)d,NF=%(fams)d,UN=%(un)d,SK=UN?'csu':'cs';var KS=['q','src','sp','ql','mn','mx','od'].concat(UN?['shp']:[]),TS=['t0','t1','t2'];function save(){try{var o={F:F,RG:$('rg').checked,IND:$('ind').checked,T:TS.map(function(i){return $(i).checked})};KS.forEach(function(i){o[i]=$(i).value});sessionStorage.setItem(SK,JSON.stringify(o))}catch(e){}}
try{var o0=JSON.parse(sessionStorage.getItem(SK)||'null');if(o0){KS.forEach(function(i){if(o0[i]!==undefined)$(i).value=o0[i]});F=o0.F||'';if(o0.T)TS.forEach(function(i,j){$(i).checked=!!o0.T[j]});$('rg').checked=!!o0.RG;$('ind').checked=o0.IND===undefined?true:!!o0.IND}}catch(e){}
var hm=(location.hash.match(/^#f=([0-9a-f]+)/)||[])[1];if(hm){F=hm;$('ind').checked=false;TS.forEach(function(i){$(i).checked=true});['q','src','sp','ql','mn','mx'].concat(UN?['shp']:[]).forEach(function(i){$(i).value=''})}
function TW(d){return ' <span class="tag tier" title="'+tr(TT[d[23]][1])+'">'+tr(TT[d[23]][0])+'</span>'+(d[26]>0?' <span class="tag ind" title="'+tr('independent structure: not contained in a larger module of the same source map')+'">'+tr('independent')+'</span>':'<a class="tag cont" href="module.html?id='+d[28]+'" title="'+tr('contained in a larger module of the same source map (overlapping window)')+'">'+tr('contained in {c}',{c:d[27]})+'</a>')+(d[24]?' <span class="tag grp" title="'+tr('same reactions drawn in more than one source map')+'">'+tr('same reactions · {n} maps',{n:d[24]})+'</span>':'')}
function UT(d){return ' <span class="tag">'+tr(d[29])+'</span>'+(d[30]?' <span class="tag">'+tr('{n} appendix reactions',{n:d[30]})+'</span>':'')}
function run(){var shp=UN?$('shp').value:'',t=$('q').value.trim().toLowerCase(),s=$('src').value,p=$('sp').value,ql=$('ql').value,mn=+$('mn').value||0,mx=+$('mx').value||1e9,od=$('od').value,r=[];
var ind=$('ind').checked,rgc=$('rg').checked,agg={sh:0,hc:0,ht:0},fams={};
D.forEach(function(d,i){if((!s||d[2]===s)&&(!p||d[3]===p)&&(!ql||d[4]===ql)&&d[5]>=mn&&d[5]<=mx&&(!t||d[11].indexOf(t)>=0)&&(!F||d[9]===F)&&(!rgc||d[24])&&(!shp||d[29]===shp)){var v=ind?d[26]:d[6];if(ind)agg.hc+=d[6]-d[26];if(!$('t'+d[23]).checked){agg.ht+=v}else if(v>0){d.i=i;r.push(d);agg.sh+=v;fams[d[9]]=1}}});
if(od==='size')r.sort(function(a,b){return b[5]-a[5]||a.i-b.i});else if(od==='name')r.sort(function(a,b){return a[1].localeCompare(b[1])||a.i-b.i});
var hid=[];if(agg.hc)hid.push(tr('{n} contained in a larger module',{n:agg.hc.toLocaleString()}));if(agg.ht)hid.push(tr('{n} by size tier',{n:agg.ht.toLocaleString()}));
$('cnt2').textContent=tr('Showing {a} of {b} module instances',{a:agg.sh.toLocaleString(),b:TOTAL.toLocaleString()})+(hid.length?' ('+tr('hidden')+': '+hid.join('; ')+')':'')+' · '+tr('counts follow the filters above; the summary cards count all modules');
$('cnt3').textContent=tr('In the current selection: {f} families, {v} layout variants (over all modules: {F} families, {V} layout variants)',{f:Object.keys(fams).length.toLocaleString(),v:r.length.toLocaleString(),F:NF.toLocaleString(),V:D.length.toLocaleString()});
cnt.textContent=tr('{a} of {b} layout variants match',{a:r.length.toLocaleString(),b:D.length.toLocaleString()})+(r.length>LIM?tr(' (first {n} shown)',{n:LIM}):'');
var prev='';out.innerHTML=r.slice(0,LIM).map(function(d){var same=d[14]===prev;prev=d[14];var rx=d[13]?esc(d[12]):'<span class="warn">'+tr('geometry only, no reaction IDs')+'</span>';
return '<tr'+(same?' class="same"':'')+'><td><a href="module.html?id='+d[21]+'">'+th(d[17])+'</a></td><td><a href="module.html?id='+d[21]+'">'+'<span class="dn">'+esc(d[1])+'</span></a>'+(UN?UT(d):TW(d))+'<div class="sub"><span class="muted" title="'+tr('module id')+' '+d[0]+'">'+d[21]+'</span> · '+'<span title="'+tr('source map id')+' '+d[22]+'">'+maptxt(d[16])+'</span> · '+esc(d[3])+'</div><div class="rx">'+rx+'</div></td><td><a href="#" class="fam" data-f="'+d[9]+'" style="--h:'+d[18]+'" title="'+tr('Click to show only this family')+'">'+d[14]+'-v'+d[15]+'</a><div class="sub">'+(d[20]>1?tr('Variant {a} of {b}',{a:d[15],b:d[20]})+' · '+tr('Same connectivity shape: {a} instances · {b} layout variants',{a:d[19],b:d[20]}):tr('Only layout')+' · '+tr('Same connectivity shape: {a} instances · only layout',{a:d[19]}))+'</div></td><td>'+esc(d[2])+'</td><td>'+esc(d[4])+'</td><td>'+d[5]+'</td><td>'+d[6]+'</td><td>'+esc(d[10])+'</td></tr>'}).join('');
save();var fc=$('fchip');if(F){fc.style.display='';$('ftxt').textContent=tr('Showing one family only')+': '+D.filter(function(x){return x[9]===F})[0][14]}else fc.style.display='none'}
out.onclick=function(e){var a=e.target.closest('a.fam');if(a){e.preventDefault();F=a.dataset.f;run();$('hits').scrollIntoView()}};
$('fclear').onclick=function(){F='';run()};document.addEventListener('langchange',run);
['q','mn','mx'].forEach(function(i){$(i).oninput=run});TS.forEach(function(i){$(i).onchange=run});$('rg').onchange=run;$('ind').onchange=run;$('tall').onclick=function(){TS.forEach(function(i){$(i).checked=true});$('ind').checked=false;run()};['src','sp','ql','od'].forEach(function(i){$(i).onchange=run});
$('reset').onclick=function(){['q','mn','mx'].forEach(function(i){$(i).value=''});['src','sp','ql'].concat(UN?['shp']:[]).forEach(function(i){$(i).value=''});$('od').value='family';F='';TS.forEach(function(i,j){$(i).checked=UN||j>0});$('rg').checked=false;$('ind').checked=!UN;run()};run()})();
</script>'''

SHAPE_FILTER='<label>Shape<select id="shp"><option value="">all</option><option value="chain">chain</option><option value="cycle">cycle</option><option value="fork">fork</option></select></label>'
UNIT_HIDDEN='<span hidden><input type="checkbox" id="t0" checked><input type="checkbox" id="t1" checked><input type="checkbox" id="t2" checked><input type="checkbox" id="ind"><input type="checkbox" id="rg"><button id="tall" type="button"></button></span>'

def make(library='windows'):
 """library 'windows' -> results/catalog (Atlas v0.3 window library); 'units' -> results/catalog_units (non-overlapping structure units)."""
 assert library in ('windows','units');U=library=='units';CUR[0]=library
 import units_library
 mods=units_library.load() if U else database.modules();families=defaultdict(list)
 for m in mods:families[m['core_family']].append(m)
 codes=assign_codes(mods,families,*(('UF','U') if U else ('F','M')));GROUPS.clear();RED.clear()
 if U:groups,gcode,redinfo,redraw={},{},{},{}  # the window library's derived indexes are keyed by window ids and must not be applied to unit ids
 else:groups,gcode=load_groups(mods,codes);redinfo,redraw=load_redundancy(mods,codes)
 GROUPS.update(groups);RED.update(redinfo)
 dest=OUT/('catalog_units' if U else 'catalog');dest.mkdir(exist_ok=True);rows=[];index=[]
 for family,instances in sorted(families.items()):
  variants={}
  for m in instances:variants.setdefault(m['variant'],m)
  for variant,m in sorted(variants.items()):
   vno=codes['v'][variant]
   count=sum(v['variant']==variant for v in instances);source_support=len({v['lineage'] for v in instances if v['variant']==variant});row={k:m[k] for k in ('id','core_family','variant','source','species','quality','size','map_id','lineage','source_sha256','license')};row.update(instances=count,source_layouts=source_support,human_review='not_required');index.append(row);rows.append((m,count,source_support,row,vno,len(instances),len(variants)))
 rows.sort(key=lambda r:(codes['f'][r[0]['core_family']],r[4]))
 pages=math.ceil(len(rows)/PER_PAGE)
 def pager(current):
  parts=[]
  for i in range(pages):parts.append(f'<span aria-current="page">{i+1}</span>' if i+1==current else f'<a href="page-{i+1}.html">{i+1}</a>')
  return '<div class="pager">'+f'<a href="index.html">Overview</a>'+''.join(parts)+'</div>'
 for i in range(pages):
  chunk=rows[i*PER_PAGE:(i+1)*PER_PAGE];cards=''.join(card(*r) for r in chunk)
  body=('<div class="note">Solid points are observed in the source layout; pale points are derived. Circles = metabolites, lines = reactions (as in the expert maps). Families ignore side attachments; chemical equivalence is a separate index.</div>'
   '<div class="filters"><label>Filter this page<input id="q" type="search" placeholder="name, species, family or reaction id"></label><label>Source<select id="src"><option value="">all</option>'+''.join(f'<option>{s}</option>' for s in SRC_COLOR)+'</select></label><span id="cnt" class="count"></span></div>'
   +pager(i+1)+f'<div class="grid">{cards}</div>'+pager(i+1))
  ui.write(dest/f'page-{i+1}.html',(f'Catalog page {i+1} of {pages}',f'图录第 {i+1} / {pages} 页') if not U else (f'Structure-unit catalog page {i+1} of {pages}',f'结构单元图录第 {i+1} / {pages} 页'),body,section='catalog',zh=ZH,script=PAGE_JS,lead=ui.bi(f'Variants {i*PER_PAGE+1}–{i*PER_PAGE+len(chunk)} of {len(rows):,}, ordered by family code (F0001 = largest family).',f'第 {i*PER_PAGE+1}–{i*PER_PAGE+len(chunk)} 个变体，共 {len(rows):,} 个，按家族编号排序（F0001 为最大家族）。'),crumbs=[('Catalog','index.html'),((f'Page {i+1}',f'第 {i+1} 页'),None)])
 with gzip.open(dest/'modules.json.gz','wt') as f:json.dump(mods,f,ensure_ascii=False)
 dump(dest/'index.json',index)
 data=[[m['id'],m['name'],m['source'],m['species'],m['quality'],m['size'],c,s,i//PER_PAGE+1,m['core_family'],m['license'],' '.join([m['name'],m['species'],m['source'],m['core_family'],m['id'],codes['f'][m['core_family']],codes['f'][m['core_family']]+'-v%d'%v,codes['m'][m['id']]]+([m['shape_type']] if U else [])+rx_ids(m)).lower(),rx_summary(m),len(m['reactions']),codes['f'][m['core_family']],v,map_line(m),thumb(m),hue(m['core_family']),fi,fv,codes['m'][m['id']],m['map_id'],tier(m['size']),(groups[m['id']]['nmaps'] if m['id'] in groups and groups[m['id']]['cross'] else 0),groups[m['id']]['code'] if m['id'] in groups else '',(redinfo[m['id']]['mx'] if m['id'] in redinfo else c),redinfo[m['id']]['cont'] if m['id'] in redinfo else '',redinfo[m['id']]['link'] if m['id'] in redinfo else '']+([m['shape_type'],len(m['appendix_reactions'])] if U else []) for i,(m,c,s,r,v,fi,fv) in enumerate(rows)]
 (dest/'search-index.js').write_text('window.CATALOG='+json.dumps(data,ensure_ascii=False,separators=(',',':'))+';')
 gallery_detail.write(dest,[(m,codes['m'][m['id']],codes['f'][m['core_family']],v,fi,fv,c,groups.get(m['id']),redinfo.get(m['id'])) for (m,c,s,r,v,fi,fv) in rows],lib=library)
 tiers=Counter(tier(m['size']) for m in mods);tier_counts=dict(total=len(mods),micro=tiers[0],small=tiers[1],structure=tiers[2])
 rs={}
 if redraw:
  mm=[m for m in mods if redraw[m['id']]['maximal']]
  rs=dict(maximal=len(mm),contained=len(mods)-len(mm),families=len({m['core_family'] for m in mm}),variants=len({m['variant'] for m in mm}),micro=sum(m['size']<=1 for m in mm),small=sum(2<=m['size']<=3 for m in mm),structure=sum(m['size']>=4 for m in mm))
 disp=dict(redundancy=rs,tiers=tier_counts,note='Display-only codes for browsing; original ids are unchanged. Families F0001.. sorted by (instance count desc, family id); variants v1.. within a family by (instance count desc, variant id); modules M00001.. sorted by (source, family code, variant number, module id).',families={c:f for f,c in codes['f'].items()},variants={codes['f'][m['core_family']]+'-v%d'%codes['v'][m['variant']]:m['variant'] for m in mods},modules={c:i for i,c in sorted(((i,c) for i,c in codes['m'].items()),key=lambda t:t[1])[::1]} if False else {c:i for i,c in sorted(((i,c) for i,c in codes['m'].items()),key=lambda t:t[1])})
 disp['groups']={c:i for i,c in sorted((c,i) for i,c in gcode.items())};disp['families']=dict(sorted(disp['families'].items()));disp['variants']=dict(sorted(disp['variants'].items()));disp['modules']=dict(sorted(((c,i) for i,c in codes['m'].items())))
 (dest/'display_ids.json').write_text(json.dumps(disp,ensure_ascii=False,separators=(',',':'),sort_keys=False))
 by_inst=Counter(m['source'] for m in mods);by_var=Counter(r['source'] for r in index);species=Counter(r['species'] for r in index);quality=Counter(m['quality'] for m in mods);licences=sorted({r['license'] for r in index})
 singles=sum(len(v)==1 for v in families.values())
 if U:
  stats='<div class="stats">'+ui.stat(f'{len(mods):,}','structure units','结构单元')+ui.stat(f'{len(families):,}','structural families')+ui.stat(f'{len(index):,}','layout variants')+ui.stat(f'{pages}','catalog pages')+ui.stat(f'{len(species)}','species')+ui.stat(f'{singles:,}','singleton families')+'</div>'
  stats+='<p class="legendline">'+ui.bi('Units never overlap, so there is no redundancy index and no size-tier filter; the table follows the filters below.','单元互不重叠，所以没有冗余索引，也没有规模分级筛选；下方表格遵循筛选条件。')+'</p>'
 else:
  stats='<div class="stats">'+ui.stat(f'{len(mods):,}','module instances')+ui.stat(f'{len(families):,}','structural families')+ui.stat(f'{len(index):,}','layout variants')+ui.stat(f'{pages}','catalog pages')+ui.stat(f'{len(species)}','species')+ui.stat(f'{singles:,}','singleton families')+'</div>'
  stats+='<p class="legendline">'+ui.bi('Summary cards count all modules; the table below follows the filters (by default: independent structures only, micro hidden).','汇总卡片统计全部模块；下方表格遵循筛选条件（默认：只显示独立结构，隐藏微型）。')+'</p>'
 breakdown=('<div class="grid"><section class="panel"><h3>Module instances by source</h3>'+bar(by_inst)+'</section><section class="panel"><h3>Layout variants by source</h3>'+bar(by_var)+'</section><section class="panel"><h3>Provenance and licence</h3><p class="meta">'+ui.bi('; '.join(f'{k} {v:,} instances' for k,v in quality.most_common())+f'. Licences: {", ".join(licences)}. Every card lists its source, licence and map id; the full record is in modules.json.gz.','; '.join(f'{k} {v:,} 个实例' for k,v in quality.most_common())+f'。许可：{", ".join(licences)}。每张卡片都列出来源、许可和来源图；完整记录见 modules.json.gz。')+'</p></section></div>')
 filters=('<h2>Search the variants</h2><div class="filters"><label>Search<input id="q" type="search" placeholder="name, species, family id or reaction id"></label>'
  '<label>Source<select id="src"><option value="">all</option>'+''.join(f'<option value="{s}">{s} ({by_var[s]:,})</option>' for s in SRC_COLOR)+'</select></label>'
  '<label>Species<select id="sp"><option value="">all</option>'+''.join(f'<option value="{html.escape(k)}">{html.escape(k)} ({v:,})</option>' for k,v in species.most_common())+'</select></label>'
  '<label>Quality<select id="ql"><option value="">all</option>'+''.join(f'<option value="{k}">{k}</option>' for k in sorted(quality))+'</select></label>'
  +(SHAPE_FILTER if U else '')+
  (UNIT_HIDDEN if U else '<fieldset class="tiers"><legend>Size tier</legend><label><input type="checkbox" id="t0" value="0"> <span title="micro: stored for provenance and statistics, not used for retrieval" data-title-zh="微型：仅为来源记录与统计保存，不参与检索">micro · 1 reaction</span></label><label><input type="checkbox" id="t1" value="1" checked> <span title="small: used only to fill gaps" data-title-zh="小型：只用于补缺">small · 2–3 reactions</span></label><label><input type="checkbox" id="t2" value="2" checked> <span title="structure: reusable local layout" data-title-zh="结构单元：可复用的局部布局">structure · ≥4 reactions</span></label><label class="rgf"><span><input type="checkbox" id="ind" checked> <span title="Independent structure = a module not contained in a larger module of the same source map; the others are overlapping windows" data-title-zh="独立结构 = 不被同一来源图中更大模块包含的模块；其余是重叠窗口">Independent structures only</span></span></label><button id="tall" type="button">Show all</button></fieldset><label class="rgf"><span><input type="checkbox" id="rg"> <span title="Modules whose exact set of reactions is also drawn in at least one other source map (lower bound)" data-title-zh="完全相同的反应集合在至少另一张来源图中也被画过（下限值）">reaction set drawn in ≥2 maps</span></span></label>')+'<label>Min size<input id="mn" type="number" min="0" style="min-width:90px;width:90px"></label><label>Max size<input id="mx" type="number" min="0" style="min-width:90px;width:90px"></label><label>Order<select id="od"><option value="family">Group by family</option><option value="size">Size (large first)</option><option value="name">Name A–Z</option></select></label><button id="reset" type="button">Reset</button></div>'
  '<p class="legendline">Same colour = same connectivity shape (family; reactions and metabolites can differ); variants = different layouts of that shape. <a href="'+ui.rel(dest/'index.html','index.html')+'#definitions">What do these terms mean?</a></p><div id="fchip" class="fchip" style="display:none"><span id="ftxt"></span><button id="fclear" type="button">Show all families</button></div><div id="cnt" class="count"></div><div id="cnt2" class="count"></div><div id="cnt3" class="count"></div><div class="tablewrap"><table><thead><tr><th>Layout</th><th>Module and reactions</th><th title="Family = modules with the same connectivity shape (reactions and metabolites can differ; only the connection pattern is the same). Variant = a different drawn layout of that shape.">Family · variant</th><th>Source</th><th>Quality</th><th>Size</th><th>Instances</th><th>Licence</th></tr></thead><tbody id="hits"></tbody></table></div>')
 if U:
  shp=Counter(m['shape_type'] for m in mods);breakdown=breakdown[:-6]+'<section class="panel"><h3>'+ui.bi('Units by shape','按形状统计的单元')+'</h3>'+bar(shp,SHAPE_COLOR)+'<p class="meta">'+ui.bi('chain = unbranched path; cycle = shortest rings (fused into ladders); fork = branch point with its arms. Appendix reactions are fragments of fewer than 3 reactions attached to a neighbouring unit.','链 = 无分支路径；环 = 最短环（可融合成梯形）；叉 = 分支点及其臂。附属反应 = 不足 3 个反应、并入相邻单元的碎片。')+'</p></section></div>'
 browse='<h2>Browse all pages</h2>'+pager(0).replace('<a href="index.html">Overview</a>','')
 files='<p>Downloads: <a href="modules.json.gz">complete JSON exchange (modules.json.gz)</a> · <a href="index.json">variant index (index.json)</a></p>'
 Wf=gallery_method.window_facts(mods if not U else database.modules());Uf=gallery_method.unit_facts(mods if U else units_library.load())
 head=gallery_method.switch(library);method=gallery_method.section(library,Wf,Uf)
 title='Layout atlas v0.3' if not U else ('Structure-unit reuse (experimental)','结构单元复用（实验性）')
 lead_units=ui.bi('Expert metabolic-map structure units (non-overlapping, 3–12 reactions) cut once from Escher community maps, WikiPathways and Reactome. Experimental: the tool uses pathway-fragment reuse by default.','从 Escher 社区地图、WikiPathways 和 Reactome 各切一次得到的专家代谢图结构单元（互不重叠，3–12 个反应）。实验性：工具默认使用通路片段复用。')
 ui.write(dest/'index.html',title,head+stats+method+breakdown+filters+browse+files,section='catalog',zh=ZH,script=INDEX_JS.replace('%(total)d',str(len(mods))).replace('%(micro)d',str(tiers[0])).replace('%(fams)d',str(len(families))).replace('%(un)d',str(int(U))),lead=lead_units if U else ui.bi('Expert metabolic-map modules from Escher community maps, WikiPathways and Reactome. Observed source layouts do not prove independent designers; results of the search link straight to the module card.','来自 Escher 社区地图、WikiPathways 和 Reactome 的专家代谢图模块。观测到的来源布局不能证明是独立设计者；搜索结果直接链接到对应模块卡片。'))
 return dict(instances=len(mods),families=len(families),variants=len(index),pages=pages)
if __name__=='__main__':
 import sys
 for lib in (sys.argv[1:] or ['windows']):print(lib,make(lib))  # python3 gallery.py [windows] [units]
