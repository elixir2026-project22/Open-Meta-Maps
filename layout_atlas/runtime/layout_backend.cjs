// JSON stdin/stdout bridge; no shell interpolation and no network services.
const fs=require('fs');const input=JSON.parse(fs.readFileSync(0,'utf8'));
(async()=>{
 if(input.backend==='elk'){
  const ELK=require('elkjs/lib/elk.bundled.js');const elk=new ELK();
  const result=await elk.layout(input.graph);process.stdout.write(JSON.stringify(result));
 }else if(input.backend==='fcose'){
  const cytoscape=require('cytoscape');cytoscape.use(require('cytoscape-fcose'));
  // fCoSE uses Math.random internally; seed it for reproducible CPU experiments.
  let seed=input.seed||20261105;Math.random=()=>{seed=(Math.imul(1664525,seed)+1013904223)>>>0;return seed/4294967296;};
  const cy=cytoscape({headless:true,styleEnabled:false,layout:{name:'preset'},elements:input.elements});
  if(input.debug_initial){const initial={};cy.nodes().forEach(n=>{initial[n.id()]=n.position();});process.stdout.write(JSON.stringify(initial));cy.destroy();return;}
  const options={name:'fcose',animate:false,randomize:false,quality:'proof',nodeDimensionsIncludeLabels:false,
   idealEdgeLength:70,nodeRepulsion:()=>4500,numIter:input.iterations||500,
   fixedNodeConstraint:input.anchors||[],packComponents:false,tile:false};
  cy.layout(options).run();let result={};cy.nodes().forEach(n=>{result[n.id()]=n.position();});
  process.stdout.write(JSON.stringify(result));cy.destroy();
 }else throw Error('unknown backend');
})().catch(e=>{process.stderr.write(String(e.stack));process.exit(1)});
