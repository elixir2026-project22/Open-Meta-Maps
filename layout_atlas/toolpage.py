"""Front end of the local Tool page (upload a GEM, run, browse skeleton and completed map) plus a canned mock job."""
import json
import ui
from core import OUT

ZH={"EC extended mode (functional layout suggestions, off by default)": "EC 扩展模式（功能性布局建议，默认关闭）", "EC extended mode adds functional layout suggestions: reactions that share an EC number but whose equations are not verified. These are suggestions, not claimed equivalences.": "EC 扩展模式会加入功能性布局建议：共享 EC 编号但方程未核实的反应。这些只是建议，并非声称等价。", "Where each reaction position comes from": "每个反应位置的来源", "Details: coverage by source": "详情：按来源的覆盖情况", "EC screening": "EC 筛选", "Library": "库", "Proposals kept": "保留的建议", "Rejected: incompatible": "因不相容被拒", "Verified equivalent": "已核实等价", "Insufficient evidence": "证据不足", "Integrity checks": "完整性检查", "Model integrity": "模型完整性", "Visualization coverage": "可视化覆盖", "Identity integrity": "身份完整性", "Integrity check failed: the job was stopped and its results should not be used.": "完整性检查未通过：任务已停止，结果不应使用。", "Information: reactions drawn in their actual flux direction (stored as reverse-only):": "说明：按实际通量方向绘制的反应（存储为仅反向）：", "Reactions": "反应", "Provenance": "来源依据", "Identity evidence": "身份证据", "Shared EC": "共享 EC", "Context compatibility": "上下文相容性", "Ambiguity": "歧义", "Structure": "结构", "Search": "搜索", "reaction id": "反应ID", "all": "全部", "{a} of {b} reactions": "{a} / {b} 个反应", "first {n} shown": "显示前 {n} 个", "Reaction": "反应", "complete equation": "完整方程", "reaction identifier": "反应标识符", "EC plus equation": "EC 加方程", "EC only": "仅 EC", "Showing one module only": "仅显示一个模块",'Layout tool':'布局工具','Tool':'工具','Run on your own model':'在你自己的模型上运行',
'Upload a genome-scale model (SBML .xml or BiGG/COBRA .json, at most 25 MB). Step 1 builds a skeleton using ONLY existing expert structures plus network bridges; step 2 completes everything else algorithmically.':'上传基因组规模代谢模型（SBML .xml 或 BiGG/COBRA .json，不超过 25 MB）。第 1 步只用已有的专家结构加网络桥接构建骨架；第 2 步用算法补全其余部分。',
'Start the tool server first: python3 project.py tool — then open http://127.0.0.1:8765/tool/index.html. This page cannot run a job when opened as a local file.':'请先启动工具服务器：python3 project.py tool，然后打开 http://127.0.0.1:8765/tool/index.html。以本地文件方式打开时，本页无法运行任务。',
'Cannot reach the tool server. Start it with: python3 project.py tool':'无法连接工具服务器。请运行：python3 project.py tool 启动它。',
'Mock mode: no server is used; a canned job is replayed.':'模拟模式：不使用服务器，回放一个预置任务。',
'1 · Model file':'1 · 模型文件','Drop a .xml or .json file here, or click to choose':'把 .xml 或 .json 文件拖到这里，或点击选择','No file chosen':'尚未选择文件','Remove':'移除',
'Only .xml (SBML) and .json (BiGG/COBRA) files are accepted.':'只接受 .xml（SBML）和 .json（BiGG/COBRA）文件。','The file is larger than 25 MB.':'文件超过 25 MB。',
'2 · Libraries':'2 · 模板库','Core (CC0 / CC BY)':'核心库（CC0 / CC BY）','EC channel':'EC 通道',


'3 · Run':'3 · 运行','Run':'运行','Cancel':'取消','Remove job':'删除任务','Progress':'进度','Status:':'状态：','Log':'日志','Job':'任务','Status':'状态','queued':'排队中','running':'运行中','done':'已完成','failed':'失败','cancelled':'已取消','pending':'等待','Seconds':'秒',
'Results':'结果','Overview':'概览','Skeleton':'骨架','Map':'地图','Modules':'模块','Downloads':'下载',
'Reactions covered by expert structures':'专家结构覆盖的反应','Complete equation':'方程完整','EC topology':'EC 拓扑','Geometry only':'仅几何','Bridge reactions':'桥接反应','Algorithm-completed':'算法补全','Total reactions':'反应总数',
'Coverage of the model reactions':'模型反应的覆盖情况','Annotation coverage':'注释覆盖率','Warnings':'警告','No warnings.':'没有警告。','Model':'模型','reactions':'反应','metabolites':'代谢物','Modules used':'使用的模块','Total time':'总用时',
'Skeleton: expert structures only, plus network bridges. Colour = source library; dashed = bridge.':'骨架：只含专家结构，加上网络桥接。颜色 = 来源库；虚线 = 桥接。','Completed map: the skeleton plus algorithm-completed reactions.':'补全后的地图：骨架加上算法补全的反应。',
'Drag to pan, wheel to zoom.':'拖动平移，滚轮缩放。','Drag to pan, wheel or pinch to zoom toward the pointer.':'拖动平移；滚轮或双指缩放（以指针为中心）。','Zoom in':'放大','Zoom out':'缩小','Minimap':'缩略图','Fit':'适应窗口','Not available yet.':'暂未提供。','Views':'视图','Draw this view':'绘制此视图','Open':'打开','not drawn':'未绘制','queued for drawing':'已排队绘制','drawn':'已绘制',
'Module':'模块','Source map':'来源图','Species':'物种','Licence':'许可','Size tier':'规模分级','Reactions covered':'覆盖的反应','Match level':'匹配级别','micro':'微型','small':'小型','structure':'结构单元',
'File':'文件','Type':'类型','Size':'大小','Download':'下载','Error':'错误','Source library':'来源库','Bridge':'桥接','no modules':'没有模块','Click a column header to sort.':'点击列标题排序。','Loading…':'加载中…',
'Uploading…':'上传中…','Waiting for the server…':'等待服务器…'}

PAGE_JS=r'''<script>
document.addEventListener('DOMContentLoaded',function(){
var $=function(i){return document.getElementById(i)},MOCK=/[?&]mock=1/.test(location.search),FILE=location.protocol==='file:';
function esc(x){return String(x).replace(/[&<>"]/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]})}
function zh(){return window.isZh&&window.isZh()}
function L(o,k){return (zh()&&o[k+'_zh'])||o[k+'_en']||o[k]||''}
function hb(n){return n>1048576?(n/1048576).toFixed(1)+' MB':n>1024?(n/1024).toFixed(0)+' KB':n+' B'}
var MAX=25*1048576,file=null,job=null,timer=null,tab='overview',sort={k:'code',d:1},mods=null,PAL={escher:'#2d8f6f',reactome:'#a8559a',wikipathways:'#3b77b8',ec:'#7a5fc0',brenda:'#b0603a'};
function url(name){return MOCK?'mock/'+name:'/api/jobs/'+job.id+'/files/'+name}
var lastErr=null;function showErr(o){lastErr=o||null;var e=$('err');if(!o){e.style.display='none';return}e.style.display='';e.innerHTML='<b>'+tr('Error')+'</b>: '+esc(typeof o==='string'?tr(o):L(o,'message'))}
/* ---- file picking ---- */
function setFile(f){showErr(null);if(!f){file=null}else if(!/\.(xml|json)$/i.test(f.name)){showErr('Only .xml (SBML) and .json (BiGG/COBRA) files are accepted.');return}else if(f.size>MAX){showErr('The file is larger than 25 MB.');return}else file=f;drawFile();$('run').disabled=!file||FILE&&!MOCK}
function drawFile(){$('finfo').innerHTML=file?'<b class="dn">'+esc(file.name)+'</b> · '+hb(file.size)+' <button type="button" id="frm" class="lnk">'+tr('Remove')+'</button>':tr('No file chosen');var b=$('frm');if(b)b.onclick=function(e){e.stopPropagation();$('fin').value='';setFile(null)}}
var dz=$('dz');dz.onclick=function(){$('fin').click()};$('fin').onchange=function(){setFile(this.files[0])};
['dragenter','dragover'].forEach(function(ev){dz.addEventListener(ev,function(e){e.preventDefault();dz.classList.add('on')})});
['dragleave','drop'].forEach(function(ev){dz.addEventListener(ev,function(e){e.preventDefault();dz.classList.remove('on')})});
dz.addEventListener('drop',function(e){var f=e.dataTransfer.files[0];if(f)setFile(f)});
/* ---- options ---- */
function opts(){return {libraries:{core:true,ec:$('o_ec').checked},views:'on_demand'}}
function notes(){var h='';if($('o_ec').checked)h+='<p class="note">'+tr('EC extended mode adds functional layout suggestions: reactions that share an EC number but whose equations are not verified. These are suggestions, not claimed equivalences.')+'</p>';$('lnote').innerHTML=h}
$('o_br').onchange=function(){if($('o_br').checked)$('o_ec').checked=true;notes()};$('o_ec').onchange=function(){if(!$('o_ec').checked&&$('o_br').checked)$('o_br').checked=false;notes()};
/* ---- run ---- */
function hint(){showErr(FILE?'Start the tool server first: python3 project.py tool — then open http://127.0.0.1:8765/tool/index.html. This page cannot run a job when opened as a local file.':'Cannot reach the tool server. Start it with: python3 project.py tool')}
$('run').onclick=function(){if(!file)return;showErr(null);job=null;mods=null;$('res').style.display='none';$('prog').style.display='';$('run').disabled=true;$('cancel').disabled=false;
 if(MOCK){mockRun();return}
 var fd=new FormData();fd.append('model',file);fd.append('options',JSON.stringify(opts()));$('slog').textContent=tr('Uploading…');
 fetch('/api/jobs',{method:'POST',body:fd}).then(function(r){return r.json().then(function(j){return [r.status,j]})}).then(function(a){if(a[0]!==201){showErr(a[1].error||'Cannot reach the tool server. Start it with: python3 project.py tool');$('run').disabled=false;return}job={id:a[1].id,status:a[1].status,steps:[],log:[],files:[],views:[]};poll()}).catch(function(){hint();$('run').disabled=false;$('prog').style.display='none'})};
$('cancel').onclick=function(){if(MOCK){clearTimeout(timer);job&&(job.status='cancelled');drawJob();$('run').disabled=false;return}if(job)fetch('/api/jobs/'+job.id,{method:'DELETE'}).then(function(){clearTimeout(timer);job.status='cancelled';drawJob();$('run').disabled=!file}).catch(hint)};
function poll(){fetch('/api/jobs/'+job.id).then(function(r){return r.json()}).then(function(j){job=j;drawJob();if(j.status==='running'||j.status==='queued')timer=setTimeout(poll,1000);else finished()}).catch(function(){hint();$('run').disabled=false})}
/* ---- mock ---- */
var MJ=null;function mockRun(){fetch('mock/job.json').then(function(r){return r.json()}).then(function(j){MJ=j;var i=0;job={id:'j_mock',status:'running',steps:j.steps.map(function(s){return {key:s.key,label_en:s.label_en,label_zh:s.label_zh,state:'pending',seconds:null}}),log:[],progress:0,warnings:[],summary:null,files:[],views:[]};
 (function tick(){if(job.status==='cancelled')return;if(i>0)Object.assign(job.steps[i-1],{state:'done',seconds:j.steps[i-1].seconds});if(i<job.steps.length){job.steps[i].state='running';job.step=job.steps[i].key;job.log.push({t:i*1.2,msg:(zh()?job.steps[i].label_zh:job.steps[i].label_en)+' …'});job.progress=i/job.steps.length;i++;drawJob();timer=setTimeout(tick,900)}else{if(/[?&]fail=1/.test(location.search)){Object.assign(job,{status:'failed',progress:1,files:[{name:'report.json',title_en:'Report',title_zh:'报告',kind:'json',bytes:600}],error:{code:'integrity',message_en:'Integrity check failed (identity_integrity).',message_zh:'完整性检查未通过（identity_integrity）。'}});job.steps[job.steps.length-1].state='failed'}else Object.assign(job,{status:'done',progress:1,summary:j.summary,files:j.files,views:j.views,warnings:j.warnings});job.log.push({t:i*1.2,msg:job.status});drawJob();finished()}})()})}
/* ---- job display ---- */
var ICON={pending:'○',running:'◐',done:'●',failed:'✕'};
function drawJob(){if(!job)return;var h='';(job.steps||[]).forEach(function(s){h+='<li class="st '+s.state+'"><span class="ic">'+(ICON[s.state]||'○')+'</span><span>'+esc(L(s,'label')||s.key)+'</span><span class="sec">'+(s.seconds!=null?s.seconds.toFixed(1)+' s':(s.state==='running'?'…':''))+'</span></li>'});
 $('steps').innerHTML=h;$('pbar').style.width=Math.round((job.progress||0)*100)+'%';$('jstat').textContent=tr(job.status||'queued');
 var lg=$('slog');lg.textContent=(job.log||[]).map(function(x){return '['+Number(x.t).toFixed(1)+'s] '+x.msg}).join('\n');lg.scrollTop=lg.scrollHeight;
 if(job.error&&lastErr!==job.error)showErr(job.error);if(job.status==='failed'||job.status==='cancelled'){$('run').disabled=!file;$('cancel').disabled=true}}
var report=null;
function loadReport(){report=null;if(!job||!(job.files||[]).some(function(f){return f.name==='report.json'})){drawInteg();return}fetch(MOCK&&/[?&]fail=1/.test(location.search)?'mock/report_fail.json':url('report.json')).then(function(r){return r.json()}).then(function(j){report=j;drawInteg()}).catch(function(){drawInteg()})}
function chk(x){if(x&&typeof x==='object'){if(x.ok!==undefined)return !!x.ok;if(x.passed!==undefined)return !!x.passed;return true}return x!==false}
function drawInteg(){var el=$('integ');if(!report||!report.integrity){el.innerHTML='';return}var I=report.integrity,ks=[['model_integrity','Model integrity'],['visualization_coverage','Visualization coverage'],['identity_integrity','Identity integrity']],bad=false;
 var h='<h3>'+tr('Integrity checks')+'</h3><ul class="checks">'+ks.map(function(k){var ok=chk(I[k[0]]);if(!ok)bad=true;return '<li class="'+(ok?'ok':'bad')+'"><b>'+(ok?'✓':'✕')+'</b> '+tr(k[1])+'</li>'}).join('')+'</ul>';
 if(bad)h='<div class="note err">'+tr('Integrity check failed: the job was stopped and its results should not be used.')+'</div>'+h;
 var ii=I.identity_integrity,dn=ii&&ii.direction_normalized_reverse_only;if(dn&&(dn.length===undefined?dn:dn.length)){var list=dn.length===undefined?String(dn):dn.join(', ');h+='<p class="meta">'+tr('Information: reactions drawn in their actual flux direction (stored as reverse-only):')+' <span class="dn">'+esc(list)+'</span></p>'}
 el.innerHTML=h}
function finished(){$('run').disabled=!file;$('cancel').disabled=true;loadReport();if(job.status==='done'){$('res').style.display='';showTab(tab)}}
/* ---- results ---- */
document.querySelectorAll('#tabs button').forEach(function(b){b.onclick=function(){showTab(b.getAttribute('data-t'))}});
function showTab(t){tab=t;document.querySelectorAll('#tabs button').forEach(function(b){b.setAttribute('aria-selected',b.getAttribute('data-t')===t?'true':'false')});var el=$('pane');if(!job||job.status!=='done'){el.innerHTML='';return}
 ({overview:overview,skeleton:function(){svgTab('skeleton.svg','Skeleton: expert structures only, plus network bridges. Colour = source library; dashed = bridge.',true)},map:mapTab,reactions:rxTab,modules:modTab,downloads:dlTab})[t]()}
var PROV=[['expert_derived','Expert-derived','专家来源','#2e7d32','Layout taken from an existing expert structure whose reaction identity was verified.','布局取自已有专家结构，且反应身份已核实。'],['topology_layout_suggestion','Topology-only suggestion','仅拓扑建议','#8e6bbf','Chosen by connection shape without any identity evidence.','仅按连接形状选择，没有任何身份证据。'],['functional_layout_suggestion','Functional (EC) suggestion','功能（EC）建议','#c9a227','Shared EC number, equation not verified: NOT a claimed equivalence.','共享 EC 编号，方程未经核实：并非声称等价。'],['bridge','Network bridge','网络桥接','#3d6fa3','Reaction added to connect expert structures in the network.','为连接专家结构而加入的网络桥接反应。'],['generated','Generated','算法生成','#9aa9b2','Placed by the layout algorithm only.','仅由布局算法放置。']];
function pl(k){var p=PROV.filter(function(x){return x[0]===k})[0];return p?(zh()?p[2]:p[1]):k}
function pcol(k){var p=PROV.filter(function(x){return x[0]===k})[0];return p?p[3]:'#888'}
function ptip(p){return zh()?p[5]:p[4]}
function overview(){var s=job.summary,c=s.covered||{},m=s.model,t=c.total||m.reactions||1,h='<p class="meta">'+tr('Model')+': <b class="dn">'+esc(m.name||m.id)+'</b> · '+m.reactions+' '+tr('reactions')+' · '+m.metabolites+' '+tr('metabolites')+' · '+tr('Total time')+': '+(s.seconds!=null?Number(s.seconds).toFixed(1)+' s':'—')+' · '+tr('Modules used')+': '+s.modules.count+'</p>';
 var pv=s.provenance;if(pv){var tot=0;PROV.forEach(function(p){tot+=pv[p[0]]||0});tot=tot||1;
  h+='<h3>'+tr('Where each reaction position comes from')+'</h3><div class="stats">'+PROV.map(function(p){return '<div class="stat" title="'+esc(ptip(p))+'" style="border-top:4px solid '+p[3]+'"><b>'+(pv[p[0]]||0).toLocaleString()+'</b><span>'+(zh()?p[2]:p[1])+'</span></div>'}).join('')+'</div>'
  +'<div class="bar big" role="img">'+PROV.map(function(p){return '<i style="width:'+(100*(pv[p[0]]||0)/tot)+'%;background:'+p[3]+'" title="'+(zh()?p[2]:p[1])+': '+(pv[p[0]]||0)+'"></i>'}).join('')+'</div><div class="legend">'+PROV.map(function(p){return '<em style="--c:'+p[3]+'" title="'+esc(ptip(p))+'">'+(zh()?p[2]:p[1])+'</em>'}).join('')+'</div>'}
 var ec=s.ec_screening;if(ec){h+='<h3>'+tr('EC screening')+'</h3><div class="tablewrap"><table><thead><tr><th>'+tr('Library')+'</th><th>'+tr('Proposals kept')+'</th><th>'+tr('Rejected: incompatible')+'</th><th>'+tr('Verified equivalent')+'</th><th>'+tr('Insufficient evidence')+'</th></tr></thead><tbody>'+Object.keys(ec).map(function(k){var e=ec[k]||{};return '<tr><td>'+esc(k)+'</td><td>'+(e.proposals_kept||0)+'</td><td>'+(e.proposals_rejected_incompatible||0)+'</td><td>'+(e.reactions_verified_equivalent||0)+'</td><td>'+(e.reactions_insufficient_evidence||0)+'</td></tr>'}).join('')+'</tbody></table></div>'}
 var cards=[['expert_module_reactions','Reactions covered by expert structures'],['complete_equation','Complete equation'],['ec_topology','EC topology'],['geometry_only','Geometry only'],['bridge_reactions','Bridge reactions'],['algorithm_completed','Algorithm-completed'],['total','Total reactions']];
 var seg=[['expert_module_reactions','#2d8f6f'],['bridge_reactions','#c58a1c'],['algorithm_completed','#8a9aa3']],det='<div class="stats">'+cards.map(function(k){return '<div class="stat"><b>'+(c[k[0]]||0).toLocaleString()+'</b><span>'+tr(k[1])+'</span></div>'}).join('')+'</div>';
 det+='<h3>'+tr('Coverage of the model reactions')+'</h3><div class="bar big" role="img">'+seg.map(function(g){var v=c[g[0]]||0;return '<i style="width:'+(100*v/(c.total||1))+'%;background:'+g[1]+'" title="'+tr(cards.filter(function(k){return k[0]===g[0]})[0][1])+': '+v+'"></i>'}).join('')+'</div><div class="legend">'+seg.map(function(g){return '<em style="--c:'+g[1]+'">'+tr(cards.filter(function(k){return k[0]===g[0]})[0][1])+' '+(c[g[0]]||0)+'</em>'}).join('')+'</div>';
 h+=pv?'<details class="det"><summary>'+tr('Details: coverage by source')+'</summary>'+det+'</details>':det;
 var a=s.annotation_coverage||{};h+='<h3>'+tr('Annotation coverage')+'</h3><div class="abars">'+['bigg','rhea','chebi','ec'].map(function(k){var v=a[k]||0;return '<div class="ab"><span>'+k.toUpperCase()+'</span><div class="bar"><i style="width:'+(100*v)+'%;background:var(--accent)"></i></div><span>'+Math.round(100*v)+'%</span></div>'}).join('')+'</div>';
 var w=job.warnings||[];h+='<h3>'+tr('Warnings')+'</h3>'+(w.length?'<ul>'+w.map(function(x){return '<li>'+esc(L(x,'message'))+'</li>'}).join('')+'</ul>':'<p class="meta">'+tr('No warnings.')+'</p>');$('pane').innerHTML=h}
function sanitize(txt){var d=new DOMParser().parseFromString(txt,'image/svg+xml'),r=d.documentElement;if(!r||r.nodeName==='parsererror'||r.getElementsByTagName('parsererror').length)return null;
 [].slice.call(r.querySelectorAll('script,foreignObject')).forEach(function(n){n.remove()});[r].concat([].slice.call(r.querySelectorAll('*'))).forEach(function(n){[].slice.call(n.attributes).forEach(function(a){if(/^on/i.test(a.name)||/^\s*javascript:/i.test(a.value))n.removeAttribute(a.name)})});return document.importNode(r,true)}
function panzoom(box,svg){if(!svg.getAttribute('viewBox')){svg.setAttribute('viewBox','0 0 '+(parseFloat(svg.getAttribute('width'))||1000)+' '+(parseFloat(svg.getAttribute('height'))||700))}svg.removeAttribute('width');svg.removeAttribute('height');svg.style.width='100%';svg.style.height='100%';
 var v0=svg.getAttribute('viewBox').split(/[\s,]+/).map(Number),vb=v0.slice();function set(){svg.setAttribute('viewBox',vb.join(' '))}
 function pt(e){var r=svg.getBoundingClientRect();return [vb[0]+(e.clientX-r.left)/r.width*vb[2],vb[1]+(e.clientY-r.top)/r.height*vb[3]]}
 function zoom(f,x,y){x=x===undefined?vb[0]+vb[2]/2:x;y=y===undefined?vb[1]+vb[3]/2:y;vb=[x-(x-vb[0])*f,y-(y-vb[1])*f,vb[2]*f,vb[3]*f];set()}
 svg.addEventListener('wheel',function(e){e.preventDefault();var p=pt(e);zoom(e.deltaY<0?.85:1.18,p[0],p[1])},{passive:false});var d=null;
 svg.addEventListener('pointerdown',function(e){d={x:e.clientX,y:e.clientY,v:vb.slice()};svg.setPointerCapture(e.pointerId)});
 svg.addEventListener('pointermove',function(e){if(!d)return;var r=svg.getBoundingClientRect();vb=[d.v[0]-(e.clientX-d.x)/r.width*d.v[2],d.v[1]-(e.clientY-d.y)/r.height*d.v[3],d.v[2],d.v[3]];set()});svg.addEventListener('pointerup',function(){d=null});
 return {fit:function(){vb=v0.slice();set()},zoom:zoom}}
function loadSvg(name,into,cb){into.innerHTML='<p class="meta">'+tr('Loading…')+'</p>';fetch(url(name)).then(function(r){if(!r.ok)throw 0;return r.text()}).then(function(t){var s=sanitize(t);if(!s)throw 0;into.innerHTML='';into.appendChild(s);cb(s)}).catch(function(){into.innerHTML='<p class="meta">'+tr('Not available yet.')+'</p>'})}
function svgTab(name,note,legend,extra){var bys=(job.summary.modules&&job.summary.modules.by_source)||{},lg='';
 if(legend){lg='<div class="legend">'+Object.keys(bys).map(function(k){return '<em style="--c:'+(PAL[k]||'#888')+'">'+esc(k)+' '+bys[k]+'</em>'}).join('')+'<em style="--c:#c58a1c">'+tr('Bridge')+'</em></div>'}
 $('pane').innerHTML='<p class="meta">'+tr(note)+' '+tr('Drag to pan, wheel to zoom.')+' <button type="button" id="pfit" class="lnk">'+tr('Fit')+'</button></p>'+lg+'<div class="viewer"><div class="vbox tall" id="pbox"></div></div>'+(extra||'');
 loadSvg(name,$('pbox'),function(s){var pz=panzoom($('pbox'),s);$('pfit').onclick=pz.fit})}
function mapTab(){var vs=job.views||[],ex='';if(vs.length){ex='<h3>'+tr('Views')+'</h3><div class="grid">'+vs.map(function(v){var done=v.state==='done';return '<div class="card"><h2 class="dn">'+esc(v.id)+'</h2><p class="meta">'+v.reactions+' '+tr('reactions')+' · '+tr(v.state==='done'?'drawn':(v.state==='queued'?'queued for drawing':'not drawn'))+'</p>'+(done?'<a href="'+url('views/'+v.id+'.svg')+'" target="_blank" rel="noopener">'+tr('Open')+'</a>':'<button type="button" class="vb" data-v="'+esc(v.id)+'"'+(v.state==='queued'?' disabled':'')+'>'+tr('Draw this view')+'</button>')+'</div>'}).join('')+'</div>'}
 svgTab('map.svg','Completed map: the skeleton plus algorithm-completed reactions.',false,ex);
 document.querySelectorAll('.vb').forEach(function(b){b.onclick=function(){drawView(b.getAttribute('data-v'))}})}
function drawView(id){var v=job.views.filter(function(x){return x.id===id})[0];v.state='queued';showTab('map');if(MOCK){setTimeout(function(){v.state='done';job.files.push({name:'views/'+id+'.svg',title_en:'View '+id,title_zh:'视图 '+id,kind:'svg',bytes:1800});if(tab==='map')showTab('map')},1500);return}
 fetch('/api/jobs/'+job.id+'/views/'+encodeURIComponent(id),{method:'POST'}).then(function(){var t=function(){fetch('/api/jobs/'+job.id).then(function(r){return r.json()}).then(function(j){job=j;if(tab==='map')showTab('map');if((j.views||[]).some(function(x){return x.state==='queued'}))setTimeout(t,1000)})};t()}).catch(hint)}
var rprov=null,rs={k:'reaction',d:1},rf='',rq='',modFilter='';
function rxTab(){var el=$('pane');el.innerHTML='<p class="meta">'+tr('Loading…')+'</p>';if(rprov)return drawRx();fetch(url('reaction_provenance.json')).then(function(r){if(!r.ok)throw 0;return r.json()}).then(function(j){rprov=j;drawRx()}).catch(function(){el.innerHTML='<p class="meta">'+tr('Not available yet.')+'</p>'})}
var EVK={complete_equation:'complete equation',reaction_identifier:'reaction identifier',ec_plus_equation:'EC plus equation',ec_only:'EC only'};
function drawRx(){if(tab!=='reactions')return;var rows=rprov.filter(function(r){return (!rf||r.provenance===rf)&&(!rq||String(r.reaction).toLowerCase().indexOf(rq)>=0)});
 function sv(r){var k=rs.k;if(k==='evidence')return r.identity_evidence?r.identity_evidence.kind:'';if(k==='structure')return r.structure?r.structure.code:'';return String(r[k]===undefined||r[k]===null?'':r[k]).toLowerCase()}
 rows.sort(function(a,b){var x=sv(a),y=sv(b);return (x<y?-1:x>y?1:0)*rs.d});var MAXR=300,shown=rows.slice(0,MAXR);
 var cols=[['reaction','Reaction'],['provenance','Provenance'],['evidence','Identity evidence'],['shared_ec','Shared EC'],['context_compatibility','Context compatibility'],['ambiguity','Ambiguity'],['structure','Structure']];
 var h='<div class="filters"><label>'+tr('Search')+'<input type="search" id="rxq" value="'+esc(rq)+'" placeholder="'+tr('reaction id')+'"></label><label>'+tr('Provenance')+'<select id="rxf"><option value="">'+tr('all')+'</option>'+PROV.map(function(p){return '<option value="'+p[0]+'"'+(rf===p[0]?' selected':'')+'>'+(zh()?p[2]:p[1])+'</option>'}).join('')+'</select></label><span class="count">'+tr('{a} of {b} reactions',{a:rows.length.toLocaleString(),b:rprov.length.toLocaleString()})+(rows.length>MAXR?' · '+tr('first {n} shown',{n:MAXR}):'')+'</span></div>';
 h+='<div class="tablewrap"><table><thead><tr>'+cols.map(function(c){return '<th><button type="button" class="lnk sortr" data-k="'+c[0]+'">'+tr(c[1])+(rs.k===c[0]?(rs.d>0?' ▲':' ▼'):'')+'</button></th>'}).join('')+'</tr></thead><tbody>'+shown.map(function(r){var e=r.identity_evidence,s=r.structure;return '<tr><td><code>'+esc(r.reaction)+'</code></td><td><span class="pdot" style="background:'+pcol(r.provenance)+'"></span>'+esc(pl(r.provenance))+'</td><td>'+(e?tr(EVK[e.kind]||e.kind)+(e.reason?'<div class="sub dn">'+esc(e.reason)+'</div>':''):'—')+'</td><td>'+(e&&e.shared_ec?esc(e.shared_ec):'—')+'</td><td>'+esc(r.context_compatibility===undefined||r.context_compatibility===null?'—':r.context_compatibility)+'</td><td>'+esc(r.ambiguity===undefined||r.ambiguity===null?'—':r.ambiguity)+'</td><td>'+(s?'<a href="#" class="tomod" data-c="'+esc(s.code)+'">'+esc(s.code)+'</a><div class="sub">'+esc(s.match_level||'')+'</div>':'—')+'</td></tr>'}).join('')+'</tbody></table></div>';
 $('pane').innerHTML=h;
 $('rxq').oninput=function(){rq=this.value.trim().toLowerCase();drawRx();var q=$('rxq');q.focus();q.setSelectionRange(q.value.length,q.value.length)};$('rxf').onchange=function(){rf=this.value;drawRx()};
 document.querySelectorAll('.sortr').forEach(function(b){b.onclick=function(){var k=b.getAttribute('data-k');rs=rs.k===k?{k:k,d:-rs.d}:{k:k,d:1};drawRx()}});
 document.querySelectorAll('.tomod').forEach(function(a){a.onclick=function(e){e.preventDefault();modFilter=a.getAttribute('data-c');showTab('modules')}})}
function modTab(){var el=$('pane');el.innerHTML='<p class="meta">'+tr('Loading…')+'</p>';if(mods)return drawMods();fetch(url('modules.json')).then(function(r){if(!r.ok)throw 0;return r.json()}).then(function(j){mods=j;drawMods()}).catch(function(){el.innerHTML='<p class="meta">'+tr('Not available yet.')+'</p>'})}
function drawMods(){if(tab!=='modules')return;var rows=mods.filter(function(m){return !modFilter||m.code===modFilter}).sort(function(a,b){var x=sv(a),y=sv(b);return (x<y?-1:x>y?1:0)*sort.d});function sv(m){var k=sort.k;if(k==='provenance')return Object.keys(m.provenance||{}).join(',');return k==='reactions'?(m.reactions||[]).length:(k==='size'?m.size:String(m[k]||'').toLowerCase())}
 var cols=[['code','Module'],['map_name','Source map'],['species','Species'],['license','Licence'],['tier','Size tier'],['reactions','Reactions covered'],['provenance','Provenance'],['match_level','Match level']];
 $('pane').innerHTML=(modFilter?'<p class="meta">'+tr('Showing one module only')+': <b>'+esc(modFilter)+'</b> <button type="button" class="lnk" id="mfc">'+tr('Show all')+'</button></p>':'')+(rows.length?'<p class="meta">'+tr('Click a column header to sort.')+'</p>':'<p class="meta">'+tr('no modules')+'</p>')+'<div class="tablewrap"><table><thead><tr>'+cols.map(function(c){return '<th><button type="button" class="lnk sortb" data-k="'+c[0]+'">'+tr(c[1])+(sort.k===c[0]?(sort.d>0?' ▲':' ▼'):'')+'</button></th>'}).join('')+'</tr></thead><tbody>'+rows.map(function(m){return '<tr><td><a href="../catalog/module.html?id='+esc(m.code)+'">'+esc(m.code)+'</a></td><td><span class="dn">'+esc(m.map_name||'')+'</span><div class="sub">'+esc(m.source||'')+'</div></td><td><span class="dn">'+esc(m.species||'')+'</span></td><td>'+esc(m.license||'')+'</td><td>'+tr(m.tier||'structure')+' · '+(m.size||'')+'</td><td title="'+esc((m.reactions||[]).join(', '))+'">'+(m.reactions||[]).length+'</td><td>'+Object.keys(m.provenance||{}).map(function(k){return '<span class="pdot" style="background:'+pcol(k)+'"></span>'+esc(pl(k))+' '+m.provenance[k]}).join('<br>')+'</td><td>'+esc(m.match_level||'')+'</td></tr>'}).join('')+'</tbody></table></div>';
 var mc=$('mfc');if(mc)mc.onclick=function(){modFilter='';drawMods()};document.querySelectorAll('.sortb').forEach(function(b){b.onclick=function(){var k=b.getAttribute('data-k');sort=sort.k===k?{k:k,d:-sort.d}:{k:k,d:1};drawMods()}})}
function dlTab(){$('pane').innerHTML='<div class="tablewrap"><table><thead><tr><th>'+tr('File')+'</th><th>'+tr('Type')+'</th><th>'+tr('Size')+'</th><th></th></tr></thead><tbody>'+(job.files||[]).map(function(f){return '<tr><td>'+esc(L(f,'title')||f.name)+'<div class="sub"><code>'+esc(f.name)+'</code></div></td><td>'+esc(f.kind)+'</td><td>'+hb(f.bytes||0)+'</td><td><a href="'+url(f.name)+'" download>'+tr('Download')+'</a></td></tr>'}).join('')+'</tbody></table></div>'}
document.addEventListener('langchange',function(){showErr(lastErr);drawInteg();if(job){drawJob();drawFile();if(job.status==='done')showTab(tab)}else drawFile()});
if(MOCK)$('mockbar').style.display='';else if(FILE)hint();
drawFile();$('run').disabled=true;
});
</script>'''

BODY='''<div class="note" id="mockbar" style="display:none">Mock mode: no server is used; a canned job is replayed.</div>
<div class="note err" id="err" style="display:none"></div>
<p class="lead2">Upload a genome-scale model (SBML .xml or BiGG/COBRA .json, at most 25 MB). Step 1 builds a skeleton using ONLY existing expert structures plus network bridges; step 2 completes everything else algorithmically.</p>
<section class="panel"><h2>1 · Model file</h2>
<div id="dz" class="dz" tabindex="0" role="button"><input type="file" id="fin" accept=".xml,.json" hidden><p>Drop a .xml or .json file here, or click to choose</p></div>
<p class="meta" id="finfo"></p></section>
<section class="panel"><h2>2 · Libraries</h2><div class="opts">
<label><input type="checkbox" id="o_core" checked disabled> <span>Core (CC0 / CC BY)</span></label>
<label><input type="checkbox" id="o_ec"> <span>EC extended mode (functional layout suggestions, off by default)</span></label>
</div><div id="lnote"></div></section>
<section class="panel"><h2>3 · Run</h2><div class="vbar2"><button type="button" id="run" class="primary">Run</button><button type="button" id="cancel" disabled>Cancel</button></div>
<div id="prog" style="display:none"><p class="meta">Status: <b id="jstat"></b></p><div class="bar big"><i id="pbar" style="width:0;background:var(--accent)"></i></div><ol class="steps" id="steps"></ol><h3>Log</h3><pre id="slog" class="slog"></pre><div id="integ"></div></div></section>
<section class="panel" id="res" style="display:none"><h2>Results</h2><div id="tabs" class="tabs" role="tablist"><button type="button" data-t="overview" aria-selected="true">Overview</button><button type="button" data-t="skeleton">Skeleton</button><button type="button" data-t="map">Map</button><button type="button" data-t="reactions">Reactions</button><button type="button" data-t="modules">Modules</button><button type="button" data-t="downloads">Downloads</button></div><div id="pane"></div></section>'''

def mock(dest):
    m=dest/'mock';(m/'views').mkdir(parents=True,exist_ok=True)
    steps=[('parse','Parse model','解析模型',0.8),('map','Map to expert library','映射到专家库',2.1),('select','Select modules','选择模块',1.4),('skeleton','Build skeleton','构建骨架',0.9),('complete','Complete layout','补全布局',6.2),('report','Write report','生成报告',0.5)]
    job=dict(id='j_mock',status='done',progress=1,steps=[dict(key=k,label_en=e,label_zh=z,state='done',seconds=s) for k,e,z,s in steps],
     summary=dict(model=dict(id='e_coli_core',name='e_coli_core (mock)',reactions=95,metabolites=72),annotation_coverage=dict(bigg=1.0,rhea=.72,chebi=.64,ec=.41),provenance=dict(expert_derived=62,topology_layout_suggestion=6,functional_layout_suggestion=3,bridge=5,generated=19),ec_screening=dict(open_core=dict(proposals_kept=3,proposals_rejected_incompatible=2,reactions_verified_equivalent=1,reactions_insufficient_evidence=4)),
      covered=dict(expert_module_reactions=70,complete_equation=68,reference_supported=66,ec_topology=12,geometry_only=2,bridge_reactions=5,algorithm_completed=20,total=95),
      modules=dict(count=11,by_source=dict(escher=8,reactome=2)),seconds=11.9),
     warnings=[dict(code='mock',message_en='This is canned mock data, not a real run.',message_zh='这是预置的模拟数据，不是真实运行结果。'),dict(code='ec_extended',message_en='EC extended mode is on: functional suggestions are not verified equivalences.',message_zh='EC 扩展模式已开启：功能性建议不是已核实的等价关系。')],
     files=[dict(name='skeleton.svg',title_en='Skeleton (SVG)',title_zh='骨架（SVG）',kind='svg',bytes=2300),dict(name='map.svg',title_en='Completed map (SVG)',title_zh='补全地图（SVG）',kind='svg',bytes=2600),dict(name='modules.json',title_en='Modules used',title_zh='使用的模块',kind='json',bytes=900),dict(name='reaction_provenance.json',title_en='Reaction provenance',title_zh='反应来源依据',kind='json',bytes=5200),dict(name='report.json',title_en='Report',title_zh='报告',kind='json',bytes=700)],
     views=[dict(id='view-1',reactions=40,state='none'),dict(id='view-2',reactions=55,state='none')],log=[])
    (m/'job.json').write_text(json.dumps(job,ensure_ascii=False,indent=1))
    mods=[dict(code='M03781',module_id='108ec78eb57ebc2a6863390f',source='escher',map_name='iJO1366.Nucleotide and histidine biosynthesis',species='Escherichia coli',license='CC0-1.0',size=1,tier='micro',reactions=['GMPR'],match_level='equation'),
          dict(code='M07702',module_id='843ef719e701a0628fcd4ad1',source='escher',map_name='iAF1260.Amino_Acid_Metabolism',species='Escherichia coli',license='CC0-1.0',size=40,tier='structure',reactions=['METS','HCYSMT2','CYSTL'],match_level='equation'),
          dict(code='M11795',module_id='',source='reactome',map_name='Selenoamino acid metabolism',species='Homo sapiens',license='CC0-1.0-data',size=8,tier='structure',reactions=['R-HSA-2408500','R-HSA-2408536'],match_level='geometry')]
    for x,p in zip(mods,[dict(expert_derived=1),dict(expert_derived=30,topology_layout_suggestion=4,bridge=3),dict(expert_derived=8,functional_layout_suggestion=3)]):x['provenance']=p
    (m/'modules.json').write_text(json.dumps(mods,ensure_ascii=False,indent=1))
    rp=[]
    kinds=[('expert_derived',dict(kind='complete_equation',shared_ec=None,reason='all participants matched'),dict(structure_id='s1',code='M03781',source='escher',match_level='equation',license='CC0-1.0',partition='core')),
           ('topology_layout_suggestion',None,dict(structure_id='s2',code='M07702',source='escher',match_level='topology',license='CC0-1.0',partition='core')),
           ('functional_layout_suggestion',dict(kind='ec_only',shared_ec='2.7.1.1',reason='shares EC 2.7.1.1; equation not verified'),dict(structure_id='s3',code='M11795',source='reactome',match_level='ec',license='CC0-1.0-data',partition='ec')),
           ('bridge',None,None),('generated',None,None)]
    for i in range(24):
        k=kinds[i%5];rp.append(dict(reaction=['PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','LDH_D','ACALD','ALCD2x'][i%12]+('' if i<12 else '_b'),provenance=k[0],identity_evidence=k[1],context_compatibility='compatible' if k[2] else None,ambiguity=None,structure=k[2]))
    (m/'reaction_provenance.json').write_text(json.dumps(rp,ensure_ascii=False,indent=1))
    ok=dict(model_integrity=True,visualization_coverage=True,identity_integrity=dict(ok=True,direction_normalized_reverse_only=['FORt']))
    bad=dict(model_integrity=True,visualization_coverage=True,identity_integrity=dict(ok=False,direction_normalized_reverse_only=[]))
    (m/'report.json').write_text(json.dumps(dict(integrity=ok),indent=1));(m/'report_fail.json').write_text(json.dumps(dict(integrity=bad),indent=1))
    def svg(extra,title):
        pts=[(80,80),(240,60),(400,110),(120,230),(300,220),(460,260)]
        s=f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 560 340"><title>{title}</title><rect width="560" height="340" fill="#fff"/>'
        for a,b,d in [(0,1,0),(1,2,0),(0,3,0),(1,4,0),(2,5,1),(3,4,0),(4,5,1)]:
            s+=f'<line x1="{pts[a][0]}" y1="{pts[a][1]}" x2="{pts[b][0]}" y2="{pts[b][1]}" stroke="{"#c58a1c" if d else "#2d8f6f"}" stroke-width="3"'+(' stroke-dasharray="8 5"' if d else '')+'/>'
        for i,(x,y) in enumerate(pts):s+=f'<circle cx="{x}" cy="{y}" r="12" fill="#cfe6df" stroke="#2d8f6f" stroke-width="2"/><text x="{x+16}" y="{y+4}" font-size="12">m{i+1}</text>'
        return s+extra+'</svg>'
    (m/'skeleton.svg').write_text(svg('','mock skeleton'))
    (m/'map.svg').write_text(svg('<line x1="120" y1="230" x2="80" y2="300" stroke="#8a9aa3" stroke-width="3"/><circle cx="80" cy="300" r="10" fill="#e5eaed" stroke="#8a9aa3" stroke-width="2"/>','mock map'))
    for v in ('view-1','view-2'):(m/'views'/f'{v}.svg').write_text(svg('',f'mock {v}'))

def make():
    dest=OUT/'tool';dest.mkdir(exist_ok=True);mock(dest)
    z=dict(ZH);z.update(ui.ZH_EXTRA)
    ui.write(dest/'index.html',('Run on your own model','在你自己的模型上运行'),BODY,section='tool',zh=ZH,script=PAGE_JS,lead='')
    return str(dest/'index.html')

if __name__=='__main__':print(make())
