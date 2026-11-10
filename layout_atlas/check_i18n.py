"""Read-only i18n audit of generated HTML: EN mode must hold no CJK; the Chinese (中文) mode must hold no runs of >=3 plain English words."""
import json,re,sys
from html.parser import HTMLParser
from pathlib import Path
import os
MIN=int(os.environ.get('I18N_MIN','3'))
ROOT=Path(__file__).resolve().parent.parent
KEEP={'open','pathway','maps','escher','wikipathways','reactome','rhea','chebi','micromap','bigg','community','cc','by','nc','sa','sbml','json','svg','gem','egnn','elk','libavoid','fcose','rmse','gpml','mnxref','metanetx','python','python3','project','py','verify','gz','html','js','e','coli','ecoli'}
VOID={'meta','link','br','img','input','hr','path','circle','line','polyline'}
class P(HTMLParser):
    def __init__(s):
        super().__init__();s.stack=[];s.items=[];s.script=''
    def handle_starttag(s,t,a):
        a=dict(a)
        if t in VOID:
            s.attr(t,a);return
        s.stack.append((t,a));s.attr(t,a)
    def attr(s,t,a):
        skip=s.skipped()
        for k in('placeholder','title','aria-label','alt'):
            if a.get(k) and not skip and t!='svg':s.items.append((a.get(k),a.get('data-'+k+'-zh') or a.get(k),True))
    def handle_endtag(s,t):
        while s.stack:
            if s.stack.pop()[0]==t:break
    def skipped(s):
        for t,a in s.stack:
            if t in('script','style','pre','code','title') or 'dn' in a.get('class','').split():
                if not(t=='title' and 0):return True
        return False
    def handle_data(s,d):
        if s.stack and s.stack[-1][0]=='script':s.script+=d;return
        if not d.strip() or s.skipped():return
        zh=None
        for t,a in reversed(s.stack):
            if 'data-zh' in a:zh=a['data-zh'];break
        s.items.append((d,zh,False))
def dictionary(script):
    m=re.search(r'var Z=(\{.*?\}),O=new WeakMap',script,flags=re.S);return json.loads(m.group(1)) if m else {}
def words(text):
    out=[];run=[];best=[]
    for tok in re.findall(r"[A-Za-z][A-Za-z0-9_.'/:-]*|[^\sA-Za-z]",text):
        if re.match(r'[A-Za-z]',tok):
            if len(tok)<=2 or re.search(r'[0-9_.:/-]',tok) or tok.lower().strip(".'") in KEEP:continue
            run.append(tok)
        elif re.match(r'[　-鿿＀-￯]',tok):
            if len(run)>len(best):best=run
            run=[]
    return max(run,best,key=len)
def audit(path):
    p=P();p.feed(path.read_text());Z=dictionary(p.script);en=[];zh=[]
    for text,alt,isattr in p.items:
        t=text.strip()
        if re.search(r'[　-鿿]',t) and not alt:en.append(t[:80])
        z=alt if (alt and not isattr) else Z.get(t,t)
        if isattr:z=alt if alt!=t else Z.get(t,t)
        if isattr and re.search(r'[　-鿿]',t):en.append(t[:80])
        w=words(z)
        if len(w)>=MIN:zh.append(' '.join(w)+'   <= '+z[:70])
    return en,zh
def pages():
    r=ROOT/'layout_atlas/results'
    yield ROOT/'index.html'
    for p in sorted(r.rglob('index.html')):yield p
    for n in (1,2,77):yield r/f'catalog/page-{n}.html'
    yield r/'catalog/module.html'
    yield r/'tool/index.html'
if __name__=='__main__' and '--dom' not in sys.argv:
    bad=0
    for path in pages():
        en,zh=audit(path);bad+=len(en)+len(zh);print(f'{path.relative_to(ROOT)}: CJK-in-EN {len(en)}, English-runs-in-中文 {len(zh)}')
        if '-v' in sys.argv:
            for x in en[:15]:print('   EN-CJK:',x)
            for x in zh[:25]:print('   ZH-EN :',x)
    print('total flags',bad);sys.exit(1 if bad else 0)

# ---------------------------------------------------------------- rendered-DOM audit (headless Chrome)
CHROME='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
SNAP='''<script>(function(){
function collect(){var out=[],w=document.createTreeWalker(document.body,NodeFilter.SHOW_TEXT),x;
 while(x=w.nextNode()){var p=x.parentElement;if(!p||p.closest('script,style,pre,code,.dn,#lang,svg,title,#__out'))continue;var t=x.nodeValue.trim();if(t)out.push(t)}
 document.querySelectorAll('[title],[aria-label],[placeholder],[alt]').forEach(function(e){if(e.closest('.dn,svg,pre,code'))return;['title','aria-label','placeholder','alt'].forEach(function(a){var v=e.getAttribute(a);if(v)out.push(v)})});return out}
function act(){%(actions)s}
setTimeout(function(){try{act()}catch(e){}setTimeout(function(){var s0={lang:document.documentElement.lang,t:collect()};document.getElementById('lang').click();setTimeout(function(){var s1={lang:document.documentElement.lang,t:collect()};document.getElementById('lang').click();setTimeout(function(){var s2={lang:document.documentElement.lang,t:collect()};var o=document.createElement('script');o.type='application/json';o.id='__out';o.textContent=JSON.stringify([s0,s1,s2]);document.body.appendChild(o)},2500)},2500)},2500)},3500)})();</script>'''
DOM_PAGES=[('home','index.html','',''),
 ('catalog-index','layout_atlas/results/catalog/index.html','',"var a=document.getElementById('tall');a&&a.click();var r=document.getElementById('rg');if(r){r.checked=true;r.dispatchEvent(new Event('change'))}"),
 ('catalog-index-all','layout_atlas/results/catalog/index.html','',"var a=document.getElementById('tall');a&&a.click()"),
 ('catalog-page','layout_atlas/results/catalog/page-1.html','',''),
 ('detail-partners','layout_atlas/results/catalog/module.html','?id=M10265',''),
 ('detail-family','layout_atlas/results/catalog/module.html','?id=M03781',''),
 ('detail-geometry-only','layout_atlas/results/catalog/module.html','?id=M15930',''),
 ('detail-maximal','layout_atlas/results/catalog/module.html','?id=M07894',''),
 ('detail-contained','layout_atlas/results/catalog/module.html','?id=M10400',''),
 ('detail-single','layout_atlas/results/catalog/module.html','?id=M11908',''),
 ('composition','layout_atlas/results/composition/index.html','',''),
 ('micromap','layout_atlas/results/micromap_extension/index.html','',''),
 ('micromap-catalog','layout_atlas/results/micromap_extension/catalog/index.html','',''),
 ('geometry','layout_atlas/results/geometry_p0/e_coli_core/index.html','',''),
 ('cases','layout_atlas/results/cases/index.html','',''),
 ('model','layout_atlas/results/models/e_coli_core/all_geometry/index.html','',''),
 ('tool','layout_atlas/results/tool/index.html','','')]
def dom_run(name,rel_path,query,actions,lang):
    import subprocess,tempfile,shutil
    scratch=Path(os.environ.get('I18N_SCRATCH','/private/tmp/claude-501/i18n_dom'));scratch.mkdir(parents=True,exist_ok=True)
    src=ROOT/rel_path;text=src.read_text()
    text=text.replace('<head>','<head>',1)
    base='<base href="file://'+str(src.parent)+'/">'
    text=text.replace('<meta charset="utf-8">','<meta charset="utf-8">'+base,1).replace('</html>',SNAP%{'actions':actions}+'</html>')
    tmp=scratch/f'{name}-{lang}.html';tmp.write_text(text)
    prof=tempfile.mkdtemp(dir=scratch)
    out=subprocess.run([CHROME,'--headless=new','--disable-gpu','--no-first-run',f'--user-data-dir={prof}',f'--lang=zh-CN','--virtual-time-budget=20000','--dump-dom',f'file://{tmp}{query}'],capture_output=True,text=True,timeout=120).stdout
    shutil.rmtree(prof,ignore_errors=True)
    m=re.search(r'<script type="application/json" id="__out">(.*?)</script>',out,flags=re.S)
    import html as H
    return json.loads(H.unescape(m.group(1))) if m else None
def dom_audit():
    bad=0
    for name,rel,query,actions in DOM_PAGES:
        res=[]
        for lang in('auto',):
            snaps=dom_run(name,rel,query,actions,lang)
            if not snaps:res.append((lang,'NO SNAPSHOT',[]));bad+=1;continue
            for s in snaps:
                items=s['t'];flags=[]
                if s['lang'].startswith('zh'):
                    for t in items:
                        w=words(t)
                        if len(w)>=MIN:flags.append('ZH-EN: '+' '.join(w)+'   <= '+t[:70])
                else:
                    flags+=['EN-CJK: '+t[:80] for t in items if re.search(r'[　-鿿]',t)]
                res.append((lang+'→'+s['lang'],len(items),flags));bad+=len(flags)
        print(f'{name}:',', '.join(f'{a}: {b} nodes, {len(c)} flags' for a,b,c in res))
        for a,b,c in res:
            for f in c[:12]:print('    ',a,f)
    print('DOM total flags',bad);return bad
if __name__=='__main__' and '--dom' in sys.argv:sys.exit(1 if dom_audit() else 0)
