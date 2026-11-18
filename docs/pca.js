(function(root){
  'use strict';
  function project(rows){
    const m=rows.length,d=rows[0].length,mean=new Float64Array(d),std=new Float64Array(d);
    for(const r of rows)for(let j=0;j<d;j++)mean[j]+=r[j];
    for(let j=0;j<d;j++)mean[j]/=m;
    for(const r of rows)for(let j=0;j<d;j++){const v=r[j]-mean[j];std[j]+=v*v;}
    for(let j=0;j<d;j++)std[j]=Math.sqrt(std[j]/m)||1;
    const A=new Float64Array(m*d),C=new Float64Array(d*d);
    for(let i=0;i<m;i++){
      const off=i*d,r=rows[i];for(let j=0;j<d;j++)A[off+j]=(r[j]-mean[j])/std[j];
      for(let j=0;j<d;j++){const v=A[off+j],base=j*d;for(let k=0;k<=j;k++)C[base+k]+=v*A[off+k];}
    }
    let total=0;
    for(let j=0;j<d;j++)for(let k=0;k<=j;k++){const v=C[j*d+k]/m;C[j*d+k]=C[k*d+j]=v;if(j===k)total+=v;}
    function top(seed){
      let v=Float64Array.from({length:d},(_,j)=>Math.sin(j*.7+seed)),lambda=0;
      for(let step=0;step<200;step++){
        const w=new Float64Array(d);let norm=0;
        for(let j=0;j<d;j++){let sum=0;for(let k=0;k<d;k++)sum+=C[j*d+k]*v[k];w[j]=sum;norm+=sum*sum;}
        norm=Math.sqrt(norm);if(norm<1e-12)return {v:new Float64Array(d),lambda:0};
        let diff=0;for(let j=0;j<d;j++){w[j]/=norm;diff+=Math.abs(w[j]-v[j]);}v=w;lambda=norm;
        if(diff<1e-7)break;
      }
      return {v,lambda};
    }
    const first=top(1),x=new Float64Array(m),y=new Float64Array(m);
    for(let j=0;j<d;j++)for(let k=0;k<d;k++)C[j*d+k]-=first.lambda*first.v[j]*first.v[k];
    const second=top(2);
    for(let i=0;i<m;i++)for(let j=0;j<d;j++){x[i]+=A[i*d+j]*first.v[j];y[i]+=A[i*d+j]*second.v[j];}
    return {x,y,ev1:first.lambda/(total||1),ev2:second.lambda/(total||1)};
  }
  const api={project};if(typeof module!=='undefined'&&module.exports)module.exports=api;root.PcmPCA=api;
})(typeof globalThis!=='undefined'?globalThis:this);
