"""Shared offline page template (inline CSS, no network) for every generated HTML index."""
import html,os,re
from core import ROOT,OUT

CSS='''
:root{color-scheme:light dark;--bg:#f4f7f8;--panel:#fff;--ink:#17323d;--mute:#5b7480;--line:#dbe4e8;--accent:#0b7285;--accent-ink:#fff;--soft:#e7f1f3;--edge:#9baab4;--rx:#b66b28;--met:#167e83;--shadow:0 1px 2px rgba(20,50,60,.06),0 4px 14px rgba(20,50,60,.05)}
@media(prefers-color-scheme:dark){:root{--bg:#10191d;--panel:#182429;--ink:#e3edf0;--mute:#90a8b2;--line:#2b3b42;--accent:#52b7c9;--accent-ink:#06242b;--soft:#1f3239;--edge:#6d818b;--rx:#e09a5a;--met:#4fc1c6;--shadow:none}}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.6 system-ui,-apple-system,"Segoe UI",sans-serif}
a{color:var(--accent);text-underline-offset:2px}
a:hover{text-decoration-thickness:2px}
.top{background:var(--panel);border-bottom:1px solid var(--line);position:sticky;top:0;z-index:5}
.top .in{max-width:1240px;margin:0 auto;padding:10px 20px;display:flex;flex-wrap:wrap;align-items:center;gap:6px 18px}
.brand{font-weight:700;color:var(--ink);text-decoration:none;letter-spacing:.01em}
.top nav{display:flex;flex-wrap:wrap;gap:6px 18px}.grp{display:flex;flex-wrap:wrap;align-items:center;gap:2px 4px}.grp+.grp{border-left:1px solid var(--line);padding-left:18px}.gl{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--mute);margin-right:4px}.stage{font-size:12px;border:1px solid var(--line);border-radius:6px;padding:0 6px;margin-right:6px;color:var(--mute)}
.top nav a{padding:4px 10px;border-radius:999px;text-decoration:none;color:var(--mute);font-size:14px}
.top nav a:hover{background:var(--soft);color:var(--ink)}
.top nav a[aria-current]{background:var(--accent);color:var(--accent-ink)}
main{max-width:1240px;margin:0 auto;padding:28px 20px 56px}
h1{font-size:clamp(24px,4vw,34px);line-height:1.2;margin:0 0 6px;letter-spacing:-.01em}
h2{font-size:20px;line-height:1.3;margin:36px 0 12px}
h3{font-size:16px;margin:0 0 6px}
.lead{color:var(--mute);max-width:78ch;margin:0 0 20px;font-size:16px}
.note{border-left:3px solid var(--accent);background:var(--soft);padding:10px 14px;border-radius:0 8px 8px 0;margin:16px 0;max-width:90ch}
.crumbs{font-size:13px;color:var(--mute);margin-bottom:12px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:18px 0 8px}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px;box-shadow:var(--shadow)}
.stat b{display:block;font-size:26px;line-height:1.15;font-variant-numeric:tabular-nums}
.stat span{color:var(--mute);font-size:13px}
.bar{display:flex;height:10px;border-radius:5px;overflow:hidden;margin:8px 0 4px;background:var(--line)}
.bar i{display:block;height:100%}
.legend{display:flex;flex-wrap:wrap;gap:4px 16px;font-size:13px;color:var(--mute)}
.legend em{font-style:normal}.legend em:before{content:"";display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px;background:var(--c)}
.panel,section.panel{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:18px 20px;box-shadow:var(--shadow);margin:16px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(290px,1fr));gap:16px;margin-top:16px}
main>article,.card,main>.grid>article{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px;box-shadow:var(--shadow);overflow-wrap:anywhere;min-width:0}
main>article{margin:18px 0}
.card h2,main article h2{font-size:15px;margin:0 0 8px;line-height:1.35}
.card svg,main article svg{display:block;width:100%;height:170px;background:var(--bg);border-radius:8px;margin-bottom:10px}
.card .meta{font-size:13px;color:var(--mute);margin:0 0 8px}
.tags{display:flex;flex-wrap:wrap;gap:5px;margin:0 0 8px}
.tag{font-size:12px;padding:1px 9px;border-radius:999px;background:var(--soft);color:var(--ink);border:1px solid var(--line)}
.tag.src-escher{border-color:#2d8f6f}.tag.src-wikipathways{border-color:#3b77b8}.tag.src-reactome{border-color:#a8559a}
.tag.lic{font-family:ui-monospace,Menlo,monospace;font-size:11px}
.kv{display:grid;grid-template-columns:auto 1fr;gap:2px 12px;font-size:13px;margin:0}
.kv dt{color:var(--mute)}.kv dd{margin:0;min-width:0;overflow-wrap:anywhere}
details{font-size:13px;margin-top:8px}
summary{cursor:pointer;color:var(--accent)}
pre{white-space:pre-wrap;word-break:break-all;font:11.5px/1.5 ui-monospace,Menlo,monospace;background:var(--bg);padding:8px;border-radius:6px;margin:6px 0 0}
code{font:.92em ui-monospace,Menlo,monospace;background:var(--soft);padding:1px 6px;border-radius:5px}
.edge{stroke:var(--edge);stroke-width:1.2}.nr{fill:var(--rx)}.nm{fill:var(--met)}
.tablewrap{overflow-x:auto;background:var(--panel);border:1px solid var(--line);border-radius:12px;box-shadow:var(--shadow);margin:14px 0}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
th,td{text-align:left;padding:9px 14px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-size:12px;text-transform:uppercase;letter-spacing:.05em;color:var(--mute);background:var(--bg);white-space:nowrap}
tr:last-child td{border-bottom:0}
tbody tr:hover td{background:var(--soft)}
main img{max-width:100%;height:auto;max-height:640px;object-fit:contain;display:block;margin:8px auto;background:#fff;border-radius:8px;border:1px solid var(--line)}
.card img{max-height:360px}
main>ul{padding-left:20px}main>ul li{margin:5px 0}
.filters{display:flex;flex-wrap:wrap;gap:10px;align-items:end;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px;box-shadow:var(--shadow);margin:16px 0}
.filters label{display:flex;flex-direction:column;font-size:12px;color:var(--mute);gap:3px}
.filters input,.filters select{font:inherit;font-size:14px;color:var(--ink);background:var(--bg);border:1px solid var(--line);border-radius:8px;padding:7px 10px;min-width:150px}
.filters input[type=search]{min-width:min(280px,100%);flex:1}
.filters button{font:inherit;font-size:14px;border:1px solid var(--line);background:var(--bg);color:var(--ink);border-radius:8px;padding:7px 14px;cursor:pointer}
.count{color:var(--mute);font-size:14px;margin:6px 0}
.pager{display:flex;flex-wrap:wrap;gap:4px;margin:14px 0}
.pager a,.pager span{min-width:34px;text-align:center;padding:4px 9px;border-radius:8px;border:1px solid var(--line);background:var(--panel);text-decoration:none;font-size:13px}
.pager span{background:var(--accent);color:var(--accent-ink);border-color:var(--accent)}
.pager a:hover{background:var(--soft)}
.hit{display:none}
footer{max-width:1240px;margin:0 auto;padding:0 20px 36px;color:var(--mute);font-size:13px}
.th{width:64px;height:64px;display:block;background:var(--bg);border-radius:6px}td .sub{font-size:12px;color:var(--mute)}td .rx,.card .rx{font-size:12px;color:var(--ink);overflow-wrap:anywhere;font-family:ui-monospace,Menlo,monospace}.card .rx{margin:0 0 8px}
.fam{text-decoration:none;cursor:pointer;display:inline-block;font:12px ui-monospace,Menlo,monospace;padding:1px 8px;border-radius:6px;white-space:nowrap;color:var(--ink);background:hsl(var(--h) 60% 50%/.18);border:1px solid hsl(var(--h) 55% 50%/.55)}
tr.same td{border-top:0}
#hits td:nth-child(2){min-width:230px}
.legendline{font-size:13px;color:var(--mute);margin:6px 0}.fchip{display:inline-flex;gap:8px;align-items:center;font-size:13px;margin:6px 0}.fchip button{font:inherit;border:1px solid var(--line);background:var(--bg);color:var(--ink);border-radius:999px;padding:1px 10px;cursor:pointer}.muted{color:var(--mute);font-size:11px;font-family:ui-monospace,Menlo,monospace}.warn{color:var(--mute);font-style:italic}dl.defs dt{font-weight:600;margin-top:8px}dl.defs dd{margin:2px 0 0;color:var(--mute)}
.toolbar{margin:0 0 12px}.toolbar a{text-decoration:none}
.viewer{background:var(--panel);border:1px solid var(--line);border-radius:12px;box-shadow:var(--shadow);overflow:hidden;margin:14px 0}
.vbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center;padding:8px 12px;border-bottom:1px solid var(--line);font-size:14px}
.vbar button{font:inherit;border:1px solid var(--line);background:var(--bg);color:var(--ink);border-radius:8px;padding:3px 12px;cursor:pointer;min-width:34px}
.vbox{height:min(70vh,640px);background:var(--bg)}.gv{width:100%;height:100%;display:block;touch-action:none;cursor:grab}
.gv .redge{fill:none;stroke:var(--ink);opacity:.6;stroke-linecap:round}.gv .redge.far{stroke-dasharray:10 7}.gv .nm.far{stroke-dasharray:4 3;fill-opacity:.1}.lsn{padding:4px 12px;margin:0;border-bottom:1px solid var(--line)}.vbar select{font:inherit;font-size:14px;color:var(--ink);background:var(--bg);border:1px solid var(--line);border-radius:8px;padding:2px 8px}
.gv .arrow{fill:var(--ink);opacity:.7}.gv .nm{fill-opacity:.3;stroke:var(--met);stroke-width:2}.gv text.rl{opacity:.85}
.gv text.lab{fill:var(--ink);paint-order:stroke;stroke:var(--bg);stroke-width:.25em;stroke-linejoin:round}.nolab .gv text.lab{display:none}
.gv g.hl circle,.gv g.hl text{filter:drop-shadow(0 0 3px var(--accent))}.gv g.hl circle{stroke:var(--accent);stroke-width:2.2}
tr.hl td{background:var(--soft)!important}
.detail-grid{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:16px}
.strip{display:flex;gap:8px;overflow-x:auto;padding-bottom:6px}.sv{display:flex;flex-direction:column;align-items:center;font-size:11px;text-decoration:none;color:var(--mute);padding:4px;border:1px solid var(--line);border-radius:8px;flex:0 0 auto}.sv[aria-current]{border-color:var(--accent);color:var(--ink);background:var(--soft)}
.pgrp{margin:6px 0 10px}.pgrp .pt{font-size:13px;color:var(--mute);margin-bottom:4px}#grpBox h3{font-size:14px;margin:12px 0 4px}#grpBox{margin-top:16px}
.pager2{display:flex;justify-content:space-between;margin-bottom:8px;font-size:14px}.svgl{display:block}.ctx{margin-bottom:6px}
@media(min-width:801px){.detail-grid.one{grid-template-columns:1fr}}
@media(max-width:800px){.detail-grid{grid-template-columns:1fr}}
.tag.tier{font-size:11px}.tag.t0{opacity:.75}.tag.ind{border-color:#2d8f6f;color:#2d8f6f}.tag.cont{border-color:var(--mute);color:var(--mute);text-decoration:none}.tag.grp{border-color:var(--accent);color:var(--accent)}label.rgf span{flex-direction:row;display:flex;gap:6px;align-items:center;font-size:13px;color:var(--ink)}label.rgf input{min-width:0;width:auto}fieldset.tiers{border:1px solid var(--line);border-radius:8px;display:flex;flex-wrap:wrap;gap:4px 12px;align-items:center;padding:4px 10px;margin:0}fieldset.tiers legend{font-size:12px;color:var(--mute);padding:0 4px}fieldset.tiers label{flex-direction:row;align-items:center;gap:4px;font-size:13px;color:var(--ink)}fieldset.tiers input{min-width:0;width:auto}
.dz{border:2px dashed var(--line);border-radius:12px;padding:28px 16px;text-align:center;cursor:pointer;background:var(--bg);color:var(--mute)}.dz.on,.dz:hover{border-color:var(--accent);color:var(--ink)}.dz p{margin:0}
.opts{display:flex;flex-wrap:wrap;gap:8px 22px}.opts label{display:flex;gap:6px;align-items:center;font-size:14px}.vbar2{display:flex;gap:10px;margin-bottom:10px}
button.primary{background:var(--accent);color:var(--accent-ink);border-color:var(--accent);font-weight:600}.vbar2 button{font:inherit;font-size:14px;border:1px solid var(--line);border-radius:8px;padding:6px 18px;cursor:pointer;background:var(--bg);color:var(--ink)}.vbar2 button:disabled{opacity:.45;cursor:not-allowed}
button.lnk{font:inherit;background:none;border:0;color:var(--accent);cursor:pointer;padding:0;text-decoration:underline}th button.lnk{text-transform:inherit;letter-spacing:inherit;color:var(--mute);font-weight:inherit;text-decoration:none;font-size:inherit}
.steps{list-style:none;padding:0;margin:10px 0}.steps .st{display:grid;grid-template-columns:24px 1fr auto;gap:8px;padding:5px 0;border-bottom:1px solid var(--line);font-size:14px;color:var(--mute)}.steps .st.running{color:var(--ink);font-weight:600}.steps .st.done{color:var(--ink)}.steps .st.failed{color:#c0392b}.steps .sec{font-variant-numeric:tabular-nums}
.slog{max-height:180px;overflow:auto}.bar.big{height:14px}.abars .ab{display:grid;grid-template-columns:60px 1fr 44px;gap:8px;align-items:center;font-size:13px;margin:4px 0}.abars .bar{margin:0}
.tabs{display:flex;flex-wrap:wrap;gap:4px;border-bottom:1px solid var(--line);margin-bottom:12px}.tabs button{font:inherit;font-size:14px;background:none;border:0;border-bottom:3px solid transparent;padding:6px 14px;cursor:pointer;color:var(--mute)}.tabs button[aria-selected=true]{color:var(--ink);border-color:var(--accent);font-weight:600}
.vbox.tall{height:min(65vh,560px)}.vbox svg{touch-action:none;cursor:grab}.note.err{border-color:#c0392b;background:rgba(192,57,43,.1)}.lead2{color:var(--mute);max-width:90ch}
.pdot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px}ul.checks{list-style:none;padding:0;margin:6px 0}ul.checks li{padding:2px 0}ul.checks .ok b{color:#2e7d32}ul.checks .bad b{color:#c0392b}details.det{margin:12px 0}details.det>summary{cursor:pointer;color:var(--accent)}
@media(max-width:640px){.top{position:static}.grp+.grp{border-left:0;padding-left:0}.top nav{flex-direction:column;gap:4px}main{padding:20px 14px 40px}.top .in{padding:8px 14px}th,td{padding:8px 10px}.grid{grid-template-columns:1fr}.filters label,.filters input,.filters select{width:100%}}
'''

GROUPS=[('',[('Map Generator','layout_atlas/results/tool/index.html','tool'),('Atlas Browser','layout_atlas/results/atlas_browser/index.html','atlas'),('About / Provenance','layout_atlas/results/about/index.html','about')])]
NAV=[x for _,g in GROUPS for x in g]

TR_JS='''<script>window.tr=function(s,v){var o=s;if(v)for(var k in v)o=o.split('{'+k+'}').join(v[k]);return o};</script>'''

def rel(path,target):
    return os.path.relpath(ROOT/target,os.path.dirname(path)).replace(os.sep,'/')

def esc(x):return html.escape(str(x))

def page(path,title,body,section='',lead='',lang='en',script='',crumbs=None,footer='',zh=None):
    """Complete HTML document; `path` is the output file, used for relative site links."""
    path=str(path)
    if isinstance(title,tuple):title=title[0]
    def link(l,h,key):return f'<a href="{rel(path,h)}"'+(' aria-current="page"' if key==section else '')+f'>{l}</a>'
    links=''.join('<div class="grp">'+(f'<span class="gl">{g}</span>' if g else '')+''.join(link(*x) for x in items)+'</div>' for g,items in GROUPS)
    top=f'<header class="top"><div class="in"><a class="brand" href="{rel(path,"index.html")}">Open Meta Maps</a><nav aria-label="Site">{links}</nav></div></header>'
    crumb=''
    if crumbs:crumb='<div class="crumbs">'+' › '.join((f'<a href="{h}">'+(bi(*l) if isinstance(l,tuple) else esc(l))+'</a>' if h else (bi(*l) if isinstance(l,tuple) else esc(l))) for l,h in crumbs)+'</div>'
    lead_html=f'<p class="lead">{lead}</p>' if lead else ''
    foot=f'<footer>{footer}</footer>' if footer else ''
    return (f'<!doctype html><html lang="{lang}"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>{esc(title if title=="Open Meta Maps" else title+" · Open Meta Maps")}</title><style>{CSS}</style>{top}<main>{crumb}<h1>{esc(title)}</h1>{lead_html}{body}</main>{foot}{script}{TR_JS}</html>')

def write(path,title,body,**kw):
    path.write_text(page(path,title,body,**kw))

def bi(en,zh=None):
    """Return escaped English label (zh ignored; kept for call-site compatibility)."""
    return esc(en)

def stat(value,label,zh=None):
    return f'<div class="stat"><b>{value}</b><span>{esc(label)}</span></div>'

def table(headers,rows,raw=False):
    """rows of cells; cells are escaped unless raw=True."""
    f=(lambda c:str(c)) if raw else esc
    return '<div class="tablewrap"><table><thead><tr>'+''.join(f'<th>{esc(h)}</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{f(c)}</td>' for c in r)+'</tr>' for r in rows)+'</tbody></table></div>'

def rewrap(path,title,section='',lead=''):
    """Re-style an existing generated page (content kept) without rerunning its experiment."""
    text=path.read_text()
    if 'class="top"' in text:
        m=re.search(r'<main>(.*)</main>',text,flags=re.S);text=m.group(1)
        text=re.sub(r'<div class="crumbs">.*?</div>','',text,count=1,flags=re.S);text=re.sub(r'<h1>.*?</h1>','',text,count=1,flags=re.S)
        l=re.match(r'\s*<p class="lead">(.*?)</p>',text,flags=re.S)
        if l:lead=lead or l.group(1);text=text[l.end():]
        path.write_text(page(path,title,text,section=section,lead=lead));return
    text=re.sub(r'<style>.*?</style>','',text,flags=re.S);text=re.sub(r'<title>.*?</title>|<!doctype html>|<meta [^>]*>|</?html[^>]*>','',text)
    text=re.sub(r'<body[^>]*>','',text);text=re.sub(r'<h1>.*?</h1>','',text,count=1,flags=re.S)
    path.write_text(page(path,title,text,section=section,lead=lead))

