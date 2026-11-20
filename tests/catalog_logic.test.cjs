const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const path=require('node:path');
const root=path.resolve(__dirname,'..');
function load(file){const context={module:{exports:{}}};vm.runInNewContext(fs.readFileSync(path.join(root,file),'utf8'),context);return context.module.exports;}
const catalog=load('docs/catalog_logic.js'),PCA=load('docs/pca.js');
const records=JSON.parse(fs.readFileSync(path.join(root,'docs/data.json')));
const refs=JSON.parse(fs.readFileSync(path.join(root,'docs/aimd_references.json')));
const passed=[];
function test(name,fn){fn();passed.push(name);}
test('Existing AIMD and imported AI remain separate and exhaustive',()=>{
  assert.equal(records.length,16944);
  assert.equal(catalog.filter(records,{origin:'AIMD'}).length,2044);
  assert.equal(catalog.filter(records,{origin:'AI'}).length,14900);
});
test('System choices cascade from the selected origin and count only that origin',()=>{
  const actual=Object.fromEntries(catalog.systemOptions(records,'AI').map(x=>[x.v,x.c]));
  assert.deepEqual(actual,{GeSbTe:4416,GeSe:4197,SbTe:2676,GeTe:833,Ge:788,Sb:669,Te:669,Se:652});
  assert.equal(Object.values(actual).reduce((a,b)=>a+b,0),14900);
  assert.ok(!catalog.systemOptions(records,'AI').some(x=>x.v==='InO'||x.v==='Si'));
  assert.ok(catalog.systemOptions(records,'AIMD').some(x=>x.v==='InO'));
  assert.equal(catalog.systemOptions(records).reduce((sum,x)=>sum+x.c,0),16944);
});
test('Dropdown filters accept multi-select arrays (OR within a field, AND across fields)',()=>{
  assert.equal(catalog.filter(records,{origin:['AIMD']}).length,2044);
  assert.equal(catalog.filter(records,{origin:['AIMD','AI']}).length,16944);
  assert.equal(catalog.filter(records,{origin:[]}).length,16944);
  const pcm=catalog.filter(records,{origin:'AIMD',functionality:'PCM'}).length,ots=catalog.filter(records,{origin:'AIMD',functionality:'OTS'}).length;
  assert.equal(catalog.filter(records,{origin:'AIMD',functionality:['PCM','OTS']}).length,pcm+ots);
});
test('Facet counts apply every other filter but ignore the field being counted',()=>{
  const counts=catalog.facetCounts(records,{origin:'AIMD',functionality:['PCM']},'functionality');
  assert.equal([...counts.values()].reduce((a,b)=>a+b,0),2044);
  assert.equal(counts.get('PCM'),catalog.filter(records,{origin:'AIMD',functionality:'PCM'}).length);
  const origins=catalog.facetCounts(records,{origin:'AIMD'},'origin');
  assert.equal(origins.get('AIMD'),2044);assert.equal(origins.get('AI'),14900);
  assert.equal(catalog.systemOptions(records,{origin:['AI'],system:['GeTe']}).reduce((a,x)=>a+x.c,0),14900);
});
test('AI subcategories and model filters preserve all five labels',()=>{
  const expected={MLIP:9263,MLIP_MACE:4416,MLIP_DeePMD:833,MLIP_NEP:4014,Diffusion:4404,'Flow Matching':1233};
  for(const [method,n] of Object.entries(expected))assert.equal(catalog.filter(records,{method}).length,n);
});
test('Error ranges are inclusive, preserve zero, and exclude missing AIMD errors',()=>{
  const sample=[{generation_origin:'AI',rdf_peak_error_pct:0,elements:[]},{generation_origin:'AI',rdf_peak_error_pct:1,elements:[]},{generation_origin:'AIMD',rdf_peak_error_pct:null,elements:[]}];
  assert.equal(catalog.filter(sample,{errorMin:0,errorMax:0}).length,1);
  assert.equal(catalog.filter(sample,{errorMin:0,errorMax:1}).length,2);
  assert.ok(catalog.filter(records,{errorMax:1}).every(r=>r.generation_origin==='AI'&&r.rdf_peak_error_pct<=1));
});
test('AI filtering adds only its actual AIMD references without changing catalog counts',()=>{
  const subset=catalog.filter(records,{method:'MLIP_DeePMD'});
  const context=catalog.withReferences(subset,refs);
  assert.equal(context.catalogCount,833);assert.equal(context.referenceCount,30);assert.equal(context.items.length,863);
  assert.ok(context.items.slice(833).every(r=>r.reference_material==='GeTe'));
  assert.equal(catalog.withReferences(catalog.filter(records,{origin:'AIMD'}),refs).referenceCount,0);
  assert.equal(catalog.withReferences(catalog.filter(records,{origin:'AI'}),refs).referenceCount,770);
});
test('The reference visibility switch preserves original AIMD catalog points and the shared map input',()=>{
  const subset=catalog.filter(records,{origin:'AI',system:'GeTe'});
  const context=catalog.withReferences(subset,refs),ids=context.items.map(x=>x.sid);
  const plain=x=>JSON.parse(JSON.stringify(x));
  assert.deepEqual(plain(catalog.mapCounts(context.items,true)),{catalogCount:833,referenceCount:30,hiddenReferenceCount:0});
  assert.deepEqual(plain(catalog.mapCounts(context.items,false)),{catalogCount:833,referenceCount:0,hiddenReferenceCount:30});
  assert.equal(context.items.filter(r=>catalog.visibleInMap(r,false)).length,833);
  assert.deepEqual(context.items.map(x=>x.sid),ids);
  assert.deepEqual(plain(catalog.mapCounts(catalog.withReferences(records,refs).items,false)),{catalogCount:16944,referenceCount:0,hiddenReferenceCount:770});
  assert.equal(catalog.mapCounts(catalog.filter(records,{origin:'AIMD'}),false).catalogCount,2044);
  assert.ok(catalog.filter(records,{origin:'AIMD'}).every(r=>catalog.visibleInMap(r,false)));
});
test('Original AIMD, AI and AIMD reference have separate styling and only references can be hidden',()=>{
  const plain=x=>JSON.parse(JSON.stringify(x));
  const aimd=records.find(r=>r.generation_origin==='AIMD'),ai=records.find(r=>r.generation_origin==='AI'),ref=refs[0];
  assert.equal(catalog.modelingGroup(aimd),'AIMD');
  assert.equal(catalog.modelingGroup(ai),'AI');
  assert.equal(catalog.modelingGroup(ref),'reference');
  const colors={AIMD:'#222222',AI:'#335599',reference:'#2a9d8f'};
  assert.deepEqual([aimd,ai,ref].map(r=>catalog.modelingColor(r,colors)),Object.values(colors));
  assert.deepEqual(plain(catalog.markerStyle({generation_origin:'AI'})),{fill:true,alpha:.5});
  assert.deepEqual(plain(catalog.markerStyle({generation_origin:'AIMD',collection_role:'reference'})),{fill:true,alpha:.8});
  assert.deepEqual(plain(catalog.markerStyle(ai,.3,.7,.2)),{fill:true,alpha:.3});
  assert.deepEqual(plain(catalog.markerStyle(aimd,.3,.7,.2)),{fill:true,alpha:.7});
  assert.deepEqual(plain(catalog.markerStyle(ref,.3,.7,.2)),{fill:true,alpha:.2});
});
test('Exact element filtering and numeric error filters compose with model filters',()=>{
  const rows=catalog.filter(records,{origin:'AI',elements:['Ge','Te'],errorMin:0,errorMax:1});
  assert.ok(rows.length>0&&rows.every(r=>r.modeling_method==='MLIP_DeePMD'));
  assert.equal(catalog.filter(records,{origin:'AI',elements:['Ge','Te','Se']}).length,0);
});
test('Joint PCA preserves rank-one covariance and handles constant vectors',()=>{
  const p=PCA.project([[0,0],[1,1],[2,2],[3,3]]);
  assert.ok(Math.abs(p.ev1-1)<1e-8&&p.ev2<1e-8);
  assert.ok(Math.abs(Array.from(p.x).reduce((a,b)=>a+b,0))<1e-8);
  const c=PCA.project([[1,1],[1,1],[1,1]]);assert.ok(Array.from(c.x).every(Number.isFinite));
});
test('Spatial bond search matches the previous all-pairs distance rule',()=>{
  const atoms=Array.from({length:120},(_,i)=>({e:i%2?'Ge':'Te',x:Math.sin(i*1.7)*10,y:Math.sin(i*2.3)*9,z:Math.cos(i*.8)*8}));
  const radius=e=>e==='Ge'?1.2:1.4,expected=[];
  for(let i=0;i<atoms.length;i++)for(let j=i+1;j<atoms.length;j++){
    const a=atoms[i],b=atoms[j],d=(a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2,cut=(radius(a.e)+radius(b.e))*1.15;
    if(d<cut*cut&&d>.64)expected.push(i+','+j);
  }
  assert.deepEqual(Array.from(catalog.bondPairs(atoms,radius),p=>p.join(',')).sort(),expected.sort());
});
module.exports={passed:passed.length,tests:passed};
if(require.main===module)console.log(JSON.stringify(module.exports,null,2));
