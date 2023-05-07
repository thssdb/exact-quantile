//import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
//
//import java.util.ArrayList;
//import java.util.List;
//
//public class AmortizedCachedForMedian {
//  private final long deltaForUnsignedCompare=1L<<63;
//  final int bufferSizeLimit,compression;
//  int bufferNum;
//  IntArrayList indexList;
//  long[] clusterMinMax,clusterSize;
//  int clusterNum;
//  long totSize;
//
//  AmortizedCachedForMedian(int compression, int bufferSizeLimit) {
//    this.compression = compression;
//    this.bufferSizeLimit = bufferSizeLimit;
//    bufferNum = 0;
//    clusterNum = 0;
//    totSize = 0;
//    int clusterNumUpperBound = (compression*2+20+bufferSizeLimit);
//    clusterMinMax = new long[clusterNumUpperBound*2];
//    clusterSize = new long[clusterNumUpperBound];
//    indexList = new IntArrayList(clusterNumUpperBound);
//  }
//  private void sortCluster(){
//    indexList.clear();
//    for(int i=0;i<clusterNum;i++)indexList.add(i);
//    indexList.sortThis((x, y) -> Long.compare(clusterMinMax[x<<1]+((clusterMinMax[x<<1|1]-clusterMinMax[x<<1])>>>1,
//        clusterMinMax[y<<1]+((clusterMinMax[y<<1|1]-clusterMinMax[y<<1])>>>1)));
//  }
//
//  private void addCluster(long size, long min,long max){
//    clusterSize[clusterNum] = size;
//    clusterMinMax[clusterNum<<1] = min;
//    clusterMinMax[clusterNum<<1|1] = max;
//    clusterNum++;
//  }
//  private long calcMid(int clusterId){return clusterMinMax[clusterId << 1]+((clusterMinMax[clusterId<<1|1]-clusterMinMax[clusterId << 1])>>>1);}
//
//  private void updateFromBuffer() {
//    sortCluster();
//    int tmpClusterNum = clusterNum;
//    clusterNum = 0;
//
//    long expectedClusterSize = (totSize + compression - 1) / compression;
//    int p1 = 0, index1=indexList.get(0), p2 = 0;
//    long min=Long.MAX_VALUE, max = Long.MIN_VALUE, size = 0, minValue,maxValue, count;
//    while (p1 < bufferNum || p2 < tmpClusterNum) {
//      if(p2==tmpClusterNum||(p1<bufferNum && (bufferValueCount[index1 << 1]/*^deltaForUnsignedCompare*/) < (calcMid(p2)/*^deltaForUnsignedCompare*/))) {
//        minValue=maxValue = bufferValueCount[index1 << 1];
//        count = bufferValueCount[index1 << 1 | 1];
//        p1++;
//        index1 = p1<bufferNum?indexList.get(p1):-1;
//      } else {
//        minValue = tmpMinMax[p2 << 1];
//        maxValue = tmpMinMax[p2<<1|1];
//        count = tmpSize[p2];
//        p2++;
//      }
//      if (size + count <= expectedClusterSize) {
//        size += count;
//        min = Math.min(minValue, min);
//        max = Math.max(maxValue, max);
//      } else {
//        if(size>0)
//          addCluster(size, min, max);
//        size = count;
//        min = minValue;
//        max = maxValue;
//      }
//    }
//    addCluster(size, min, max);
//
//    bufferNum = 0;
//
////    System.out.print("\t\t\t");
////    for(int i=0;i<clusterNum;i++)System.out.print(clusterMinMax[i<<1]+"---"+clusterMinMax[i<<1|1]+":"+clusterSize[i]+"\t");
////    System.out.println();
//  }
//
//  public void add(final long value,final long count){
////    System.out.println("\t\t add: "+value+"  "+count);
//    totSize += count;
//    bufferValueCount[bufferNum<<1]=value^deltaForUnsignedCompare;
//    bufferValueCount[bufferNum<<1|1]=count;
//    bufferNum++;
//    if(bufferNum == bufferSizeLimit)
//      updateFromBuffer();
//  }
//  public List<Long> findResultRange(long K1, long K2){
//    System.out.println("[DEBUG] findResultRange   cluster:");
////    System.out.print("\t\t\t");
////    for(int i=0;i<clusterNum;i++)System.out.print((clusterMinMax[i<<1]^deltaForUnsignedCompare)+"~"+(clusterMinMax[i<<1|1]^deltaForUnsignedCompare)+":"+clusterSize[i]+"\t");
////    System.out.println();
////    System.out.print("\t\t\t");
////    for(int i=0;i<clusterNum;i++)System.out.print((clusterMinMax[i<<1|1]-clusterMinMax[i<<1])+","+clusterSize[i]+"\t");
////    System.out.println();
//    double totInterval = 0,MMPSIZE=0,FAKEMIN = 1e80;
//    for(int i=0;i<clusterNum;i++){
//      double interval=clusterMinMax[i<<1|1]-clusterMinMax[i<<1];
//      totInterval+=interval;
//      if(interval<=2333.0)
//        MMPSIZE+=clusterSize[i];
//      else
//        if(Math.abs(clusterMinMax[i<<1]^deltaForUnsignedCompare)<=2333.0)
//          FAKEMIN = Math.min(FAKEMIN, clusterMinMax[i<<1]^deltaForUnsignedCompare);
//    }
//    System.out.println("[avg interval]:"+ totInterval /clusterNum + "\t\t precise_percent:"+MMPSIZE/totSize+"   FAKE_MIN:"+FAKEMIN);
//    List<Long> result = new ArrayList<>(7);
//    if(bufferNum>0)
//      updateFromBuffer();
//    long sum=0;
//    for(int i=0;i<clusterNum;i++){
//      if(sum+clusterSize[i]>=K1&&result.size()==0) {
//        result.add(sum);
//        result.add(clusterMinMax[i<<1]^deltaForUnsignedCompare);
//        result.add(clusterMinMax[i<<1|1]^deltaForUnsignedCompare);
//      }
//      if(sum+clusterSize[i]>=K2&&result.size()==3) {
//        result.add(sum);
//        result.add(clusterMinMax[i<<1]^deltaForUnsignedCompare);
//        result.add(clusterMinMax[i<<1|1]^deltaForUnsignedCompare);
//        return result;
//      }
//      sum+=clusterSize[i];
//    }
//    return result;
//  }
//
//
//}