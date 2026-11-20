(function(root){
  'use strict';
  const origin=r=>r.generation_origin||'AIMD';
  const error=r=>typeof r.rdf_peak_error_pct==='number'&&Number.isFinite(r.rdf_peak_error_pct)?r.rdf_peak_error_pct:null;
  // 下拉筛选值：''/null/[] 表示不限；字符串为单选；数组/Set 为多选（同一维度内取并集，维度之间取交集）
  function pick(sel){
    if(sel==null||sel==='')return null;
    if(typeof sel==='string')return new Set([sel]);
    const s=new Set(sel);return s.size?s:null;
  }
  const FACETS={origin:r=>origin(r),system:r=>r.system,phase:r=>r.phase,temperature:r=>r.temperature,
    functionality:r=>r.functionality||'Others',contributor:r=>r.contributor||'none'};
  function filter(records,f={}){
    const q=String(f.q||'').trim().toLowerCase(),elements=f.elements||[];
    const sets=Object.keys(FACETS).map(k=>[FACETS[k],pick(f[k])]).filter(x=>x[1]);
    return records.filter(r=>{
      for(const [get,set] of sets)if(!set.has(get(r)))return false;
      if(f.method&&r.generation_family!==f.method&&r.modeling_method!==f.method)return false;
      if(elements.length&&(r.elements.length!==elements.length||!elements.every(e=>r.elements.includes(e))))return false;
      if(f.errorMin!=null||f.errorMax!=null){const value=error(r);
        if(value===null||(f.errorMin!=null&&value<f.errorMin)||(f.errorMax!=null&&value>f.errorMax))return false;}
      return !q||[r.formula,r.system,r.path,r.modeling_method,r.sid,r.source_file].some(x=>String(x||'').toLowerCase().includes(q));
    });
  }
  // 分面计数：某一维度的选项计数 = 应用「除该维度外」所有筛选后的结构数（同维度多选不互相清零）
  function facetCounts(records,f,key){
    const counts=new Map(),get=FACETS[key];
    filter(records,{...f,[key]:null}).forEach(r=>{const v=get(r);counts.set(v,(counts.get(v)||0)+1);});
    return counts;
  }
  function withReferences(items,references){
    const ids=new Set(),materials=new Set();
    items.filter(r=>origin(r)==='AI').forEach(r=>{(r.reference_ids||[]).forEach(id=>ids.add(id));if(r.reference_material)materials.add(r.reference_material);});
    const existing=new Set(items.map(r=>r.sid));
    const context=references.filter(r=>!existing.has(r.sid)&&(ids.has(r.reference_id)||(!ids.size&&materials.has(r.reference_material))));
    return {items:items.concat(context),catalogCount:items.length,referenceCount:context.length};
  }
  // 第二个参数可为生成方式（兼容旧调用）或完整筛选对象；计数忽略体系本身的选择
  function systemOptions(records,f={}){
    if(typeof f==='string'||Array.isArray(f))f={origin:f};
    const systems=new Map();
    filter(records,{...f,system:null}).forEach(r=>{
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
  const api={origin,error,filter,facetCounts,withReferences,systemOptions,modelingGroup,modelingColor,visibleInMap,mapCounts,markerStyle,bondPairs};
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
  root.PcmCatalog=api;
})(typeof globalThis!=='undefined'?globalThis:this);
