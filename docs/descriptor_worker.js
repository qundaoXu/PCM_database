importScripts('pca.js');
self.onmessage=event=>{
  try{
    const {matrix,dim}=event.data,rows=[];
    for(let off=0;off<matrix.length;off+=dim)rows.push(matrix.subarray(off,off+dim));
    const result=PcmPCA.project(rows);
    self.postMessage(result,[result.x.buffer,result.y.buffer]);
  }catch(error){self.postMessage({error:String(error)});}
};
