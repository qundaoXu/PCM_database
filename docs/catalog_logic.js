(function(root){
  'use strict';
  const origin=r=>r.generation_origin||'AIMD';
  const error=r=>typeof r.rdf_peak_error_pct==='number'&&Number.isFinite(r.rdf_peak_error_pct)?r.rdf_peak_error_pct:null;
  function filter(records,f={}){
    const q=String(f.q||'').trim().toLowerCase(),elements=f.elements||[];
    return records.filter(r=>{
      if(f.origin&&origin(r)!==f.origin)return false;
      if(f.method&&r.generation_family!==f.method&&r.modeling_method!==f.method)return false;
      if(f.system&&r.system!==f.system)return false;
      if(f.phase&&r.phase!==f.phase)return false;
      if(f.temperature&&r.temperature!==f.temperature)return false;
      if(f.functionality&&(r.functionality||'Others')!==f.functionality)return false;
      if(f.contributor&&(r.contributor||'none')!==f.contributor)return false;
      if(elements.length&&(r.elements.length!==elements.length||!elements.every(e=>r.elements.includes(e))))return false;
      if(f.errorMin!=null||f.errorMax!=null){const value=error(r);
        if(value===null||(f.errorMin!=null&&value<f.errorMin)||(f.errorMax!=null&&value>f.errorMax))return false;}
      return !q||[r.formula,r.system,r.path,r.modeling_method,r.sid,r.source_file].some(x=>String(x||'').toLowerCase().includes(q));
    });
  }
  function withReferences(items,references){
    const ids=new Set(),materials=new Set();
    items.filter(r=>origin(r)==='AI').forEach(r=>{(r.reference_ids||[]).forEach(id=>ids.add(id));if(r.reference_material)materials.add(r.reference_material);});
    const existing=new Set(items.map(r=>r.sid));
    const context=references.filter(r=>!existing.has(r.sid)&&(ids.has(r.reference_id)||(!ids.size&&materials.has(r.reference_material))));
    return {items:items.concat(context),catalogCount:items.length,referenceCount:context.length};
  }
  function systemOptions(records,selectedOrigin=''){
    const systems=new Map();
    filter(records,{origin:selectedOrigin}).forEach(r=>{
      if(!systems.has(r.system))systems.set(r.system,{v:r.system,n:r.elements.length,c:0});
      systems.get(r.system).c++;
    });
    return [...systems.values()].sort((a,b)=>b.c-a.c||a.v.localeCompare(b.v));
  }
  function modelingGroup(record){
    return record.collection_role==='reference'?'reference':origin(record)==='AI'?'AI':'AIMD';
  }
  function modelingColor(record,colors){return colors[modelingGroup(record)];}
  function visibleInMap(record,showReference=true){return showReference||modelingGroup(record)!=='reference';}
  function mapCounts(items,showReference=true){
    let catalogCount=0,referenceCount=0,hiddenReferenceCount=0;
    items.forEach(r=>{if(!visibleInMap(r,showReference))hiddenReferenceCount++;
      else if(r.collection_role==='reference')referenceCount++;else catalogCount++;});
    return {catalogCount,referenceCount,hiddenReferenceCount};
  }
  function markerStyle(record,aiAlpha=0.5,aimdAlpha=0.8,referenceAlpha=0.8){
    return {fill:true,alpha:{AI:aiAlpha,AIMD:aimdAlpha,reference:referenceAlpha}[modelingGroup(record)]};
  }
  function bondPairs(atoms,radius,factor=1.15){
    const width=2*Math.max(...atoms.map(a=>radius(a.e)))*factor,buckets=new Map(),cells=[];
    atoms.forEach((a,i)=>{const cell=[Math.floor(a.x/width),Math.floor(a.y/width),Math.floor(a.z/width)],key=cell.join(',');
      cells.push(cell);if(!buckets.has(key))buckets.set(key,[]);buckets.get(key).push(i);});
    const pairs=[];
    atoms.forEach((a,i)=>{const c=cells[i];
      for(let dx=-1;dx<=1;dx++)for(let dy=-1;dy<=1;dy++)for(let dz=-1;dz<=1;dz++){
        const neighbors=buckets.get([c[0]+dx,c[1]+dy,c[2]+dz].join(','));if(!neighbors)continue;
        neighbors.forEach(j=>{if(j<=i)return;const b=atoms[j],x=a.x-b.x,y=a.y-b.y,z=a.z-b.z,d2=x*x+y*y+z*z,
          cut=(radius(a.e)+radius(b.e))*factor;if(d2<cut*cut&&d2>0.64)pairs.push([i,j]);});
      }
    });return pairs;
  }
  const api={origin,error,filter,withReferences,systemOptions,modelingGroup,modelingColor,visibleInMap,mapCounts,markerStyle,bondPairs};
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
  root.PcmCatalog=api;
})(typeof globalThis!=='undefined'?globalThis:this);
