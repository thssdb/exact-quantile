public class TDigestForPageStatCalculator extends TDigestForPageStat{
  int N;
  long[] clusterMidMinMaxSize;
  int compression,maxClusterNum,clusterNum;
  double bufferRate;
  int compactTimes=0;



  public TDigestForPageStatCalculator(int maxMemoryByte,double bufferRate) {
    this.maxClusterNum = maxMemoryByte/(4*Long.BYTES);
    this.compression = (int)(maxClusterNum*(1-bufferRate))/2;
    this.bufferRate = bufferRate;
    this.N = 0;
    this.clusterMidMinMaxSize = new long[maxClusterNum*4];
    this.clusterNum=0;
  }

  public void show(){
    long tmpN=0;
    for(int i=0;i<clusterNum<<2;i+=4) {
      tmpN+=clusterMidMinMaxSize[i|3];
    }
    System.out.println("\t\t[TDigestCalculator]:\tcompression:"+compression+
        "\tmaxClusterNum:"+maxClusterNum+"\tclusterNum:"+clusterNum);
    System.out.println("\t\tN:"+N+"\ttmpN:"+tmpN);
  }
  public void showCluster(){
    sortCluster(0,clusterNum-1);
    for(int i=0;i<clusterNum<<2;i+=4) {
      System.out.println("\t\t[cluster]:\t" + '[' + clusterMidMinMaxSize[i | 1] +
          "," + clusterMidMinMaxSize[i | 2] + "]:\t" + clusterMidMinMaxSize[i | 3]);
    }
  }
  public long getN(){return N;}
  public int getCompactTimes(){return compactTimes;}


  public void addCluster(long L,long R,long size){
    N+=size;
    if(clusterNum>0&&L==clusterMidMinMaxSize[(clusterNum-1)<<2|1]&&R==clusterMidMinMaxSize[(clusterNum-1)<<2|2]){
      clusterMidMinMaxSize[(clusterNum-1)<<2|3]+=size;
      return;
    }
    if(clusterNum+1==maxClusterNum)
      compact();
    int index=(clusterNum++)<<2;
    long mid =L+((R-L)>>>1);
    clusterMidMinMaxSize[index]=mid;
    clusterMidMinMaxSize[index|1]=L;
    clusterMidMinMaxSize[index|2]=R;
    clusterMidMinMaxSize[index|3]=size;
  }
  private void reAddCluster(long L,long R,long size){
    int index=(clusterNum++)<<2;
    long mid =L+((R-L)>>>1);
    clusterMidMinMaxSize[index]=mid;
    clusterMidMinMaxSize[index|1]=L;
    clusterMidMinMaxSize[index|2]=R;
    clusterMidMinMaxSize[index|3]=size;
  }

  public void updateFromPageStat(TDigestForPageStatReader reader){
    long L,R,size;
    while(reader.hasNextCluster()){
      L=reader.readNext();
      R=reader.readNext();
      size=reader.readNext();
      addCluster(L,R,size);
    }
  }
  public void updateFromData(long num){
    addCluster(num,num,1);
  }

  private void sortCluster(int l,int r) {
    if (l >= r) return;
    if(l+6>=r){
      for(int i=0;i<r-l;i++){
        boolean flag=true;
        for(int j=l;j<r;j++)
          if(clusterMidMinMaxSize[j<<2]>clusterMidMinMaxSize[(j+1)<<2]){
            swapCluster(j,j+1);
            flag=false;
          }
        if(flag)break;
      }
      return;
    }
    swapCluster(l, (l + r + 1) >>> 1);
    long pivotV = clusterMidMinMaxSize[l << 2];
    int lt = l, gt = r + 1, i = l + 1;
    while (i < gt) {
      if (clusterMidMinMaxSize[i << 2] < pivotV) {
        swapCluster(i, lt + 1);
        i++;
        lt++;
      }else if(clusterMidMinMaxSize[i<<2]>pivotV){
        swapCluster(i,gt-1);
        gt--;
      }else i++;
    }
    swapCluster(l,lt);
    sortCluster(l,lt-1);
    sortCluster(gt,r);
  }
  private void swapCluster(int id1,int id2){
    long tmp;
    id1<<=2;id2<<=2;
    for(int i=0;i<4;i++) {
      tmp = clusterMidMinMaxSize[id1];
      clusterMidMinMaxSize[id1] = clusterMidMinMaxSize[id2];
      clusterMidMinMaxSize[id2] = tmp;
      id1++;
      id2++;
    }
  }

  public void compact(){
    compactTimes++;
    sortCluster(0,clusterNum-1);
    long expectedClusterSize = (N + compression - 1) / compression;
//    System.out.println("\t\t\t[compact]\texpSize:"+expectedClusterSize+"\tN:"+N);
    long cntMin = Long.MAX_VALUE,cntMax=Long.MIN_VALUE,cntSize=0;
    long v1,v2,v3;
    int oldClusterNum = clusterNum;
    clusterNum = 0;
    for(int i=0;i<oldClusterNum;i++){
      v1=clusterMidMinMaxSize[i<<2|1];
      v2=clusterMidMinMaxSize[i<<2|2];
      v3=clusterMidMinMaxSize[i<<2|3];
//      System.out.println("\t\t"+v1+" "+v2+":"+v3);
      if(cntSize==0||cntSize+v3<=expectedClusterSize){
        cntMin=Math.min(cntMin,v1);
        cntMax=Math.max(cntMax,v2);
        cntSize+=v3;
      }else{
        reAddCluster(cntMin,cntMax,cntSize);
        cntMin=v1;
        cntMax=v2;
        cntSize=v3;
      }
    }
    reAddCluster(cntMin,cntMax,cntSize);
  }

  public long possibleSizeLEValue(long V){
    long ans=0;
    for(int i=0;i<clusterNum<<2;i+=4) {
      if (clusterMidMinMaxSize[i | 2] <= V)
        ans += clusterMidMinMaxSize[i | 3];
      else if (clusterMidMinMaxSize[i | 1] <= V && V < clusterMidMinMaxSize[i | 2])
        ans += clusterMidMinMaxSize[i | 3] - 1;
    }
    return ans;
  }
  public long possibleSizeGEValue(long V){
    long ans=0;
    for(int i=0;i<clusterNum<<2;i+=4) {
      if (V<=clusterMidMinMaxSize[i | 1])
        ans += clusterMidMinMaxSize[i | 3];
      else if (clusterMidMinMaxSize[i | 1] < V && V <= clusterMidMinMaxSize[i | 2])
        ans += clusterMidMinMaxSize[i | 3] - 1;
    }
    return ans;
  }

  public long minValueWithRank(long K){
    long L=Long.MIN_VALUE,R=Long.MAX_VALUE,mid;
    while(L<R) {
      mid = L + ((R - L) >>> 1);
      if (possibleSizeLEValue(mid) < K) L = mid+1;
      else R = mid;
    }
    return L;
  }
  public long maxValueWithRank(long K){
    K = N-K+1;
    long L=Long.MIN_VALUE,R=Long.MAX_VALUE,mid;
    while(L<R) {
      mid = L + ((R - L) >>> 1);
      if(mid==L)mid++;
      if (possibleSizeGEValue(mid) < K) R = mid-1;
      else L = mid;
    }
    return L;
  }

  public void reset(){
    clusterNum = 0;
    N = 0;
    compactTimes=0;
  }

}