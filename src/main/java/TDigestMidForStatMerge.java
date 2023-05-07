// inspired by t-Digest by Ted Dunning. See https://github.com/tdunning/t-digest
// This is a simple implementation with radix sort and K0.
// Clusters are NOT strictly in order.

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.ArrayList;
import java.util.List;

public class TDigestMidForStatMerge {
  public final long deltaForUnsigned;
  public final int compression;
  public long[] clusterMinMax, clusterSize, tmpMinMax, tmpSize;
  public int clusterNum, clusterNumLimit;
  public int maxSeriByte, maxMemByte;
  public long totN;
  boolean sorted = true;

  public TDigestMidForStatMerge(int maxMemByte, int maxSeriByte) {
    this.maxMemByte = maxMemByte;
    this.maxSeriByte = maxSeriByte;
    this.clusterNumLimit = maxMemByte / 3 / 8;
    this.compression = this.clusterNumLimit / 2 * 2 / 4;
    clusterNum = 0;
    totN = 0;
    clusterMinMax = new long[clusterNumLimit * 2];
    clusterSize = new long[clusterNumLimit];
    deltaForUnsigned = 1L << (64 - 1);
  }

  private void addCluster(long size, long min, long max) {
    clusterSize[clusterNum] = size;
    clusterMinMax[clusterNum << 1] = min;
    clusterMinMax[clusterNum << 1 | 1] = max;
    clusterNum++;
    if (clusterNum == clusterNumLimit) compaction(this.compression);
  }

  private long calcTmpMid(int clusterId) {
    return tmpMinMax[clusterId << 1] + ((tmpMinMax[clusterId << 1 | 1] - tmpMinMax[clusterId << 1]) >>> 1);
  }
  private long calcMid(int clusterId) {
    return clusterMinMax[clusterId << 1] +
        ((clusterMinMax[clusterId << 1 | 1] - clusterMinMax[clusterId << 1]) >>> 1);
  }

  public void sortCluster() {
    tmpMinMax = new long[clusterNum * 2];
    tmpSize = new long[clusterNum];
    IntArrayList indexList = new IntArrayList(clusterNum);
    for (int i = 0; i < clusterNum; i++) indexList.add(i);
    indexList.sort((x, y) -> (Long.compare(calcMid(x), calcMid(y))));

    for (int i = 0; i < clusterNum; i++) {
      int index = indexList.getInt(i);
      tmpSize[i] = clusterSize[index];
      tmpMinMax[i << 1] = clusterMinMax[index << 1];
      tmpMinMax[i << 1 | 1] = clusterMinMax[index << 1 | 1];
    }
    clusterSize = tmpSize;
    clusterMinMax = tmpMinMax;
    tmpMinMax = tmpSize = null;
    sorted = true;
  }

  private void compaction(int compression) {
    tmpMinMax = new long[clusterNum * 2];
    tmpSize = new long[clusterNum];
    System.arraycopy(clusterMinMax, 0, tmpMinMax, 0, clusterNum * 2);
    System.arraycopy(clusterSize, 0, tmpSize, 0, clusterNum);
    int tmpClusterNum = clusterNum;

    IntArrayList indexList = new IntArrayList(clusterNum);
    for (int i = 0; i < clusterNum; i++) indexList.add(i);
    indexList.sort((x, y) -> (Long.compare(calcMid(x), calcMid(y))));
    clusterNum = 0;

    long expectedClusterSize = (totN + compression - 1) / compression;
    int index;
    long min = Long.MAX_VALUE, max = Long.MIN_VALUE, size = 0;
    for (int i = 0; i < tmpClusterNum; i++) {
      index = indexList.getInt(i);
      if(compression!=this.compression&&tmpClusterNum-i+clusterNumLimit<=maxSeriByte/3/8){
        for(;i<tmpClusterNum;i++){
          index = indexList.getInt(i);
          addCluster(tmpSize[index],tmpMinMax[index << 1],tmpMinMax[index << 1 | 1]);
        }
        break;
      }
      if (size + tmpSize[index] <= expectedClusterSize) {
        size += tmpSize[index];
        min = Math.min(min, tmpMinMax[index << 1]);
        max = Math.max(max, tmpMinMax[index << 1 | 1]);
      } else {
        if (size > 0)
          addCluster(size, min, max);
        size = tmpSize[index];
        min = tmpMinMax[index << 1];
        max = tmpMinMax[index << 1 | 1];
      }
    }
    if(size>0)
    addCluster(size, min, max);
    sorted = true;

    tmpMinMax = tmpSize = null;

//    System.out.println("\t\t expectedClusterSize:"+expectedClusterSize + " totN:"+totN+" oldClusterNum:"+tmpClusterNum);
  }

//  public void add(final long value,final long count){
////    System.out.println("\t\t add: "+(value)+"  "+count);
//    totN += count;
////    bufferValueCount[bufferNum<<1]=value;
////    bufferValueCount[bufferNum<<1|1]=count;
//    bufferNum++;
//    if(bufferNum == bufferSizeLimit)
//      updateFromBuffer();
//  }

  public void update(long value) {
//    System.out.println("\t\t add: "+(value));
    totN++;
//    value ^= deltaForUnsigned;
    sorted = false;
    addCluster(1, value, value);
  }

  public void merge(TDigestMidForStatMerge another) {
//    System.out.println("\t\t add: "+(value));
    totN += another.totN;
    for (int i = 0; i < another.clusterNum; i++)
      addCluster(another.clusterSize[i], another.clusterMinMax[i << 1], another.clusterMinMax[i << 1 | 1]);
  }
  public void merge(List<TDigestMidForStatMerge> anotherList) {
    for(TDigestMidForStatMerge another:anotherList)
      merge(another);
  }

  private long possibleSizeLEValue(final long V) {
    long size = 0;
    for (int i = 0; i < clusterNum; i++)
      if (clusterMinMax[i << 1 | 1] <= V)
        size += clusterSize[i];
      else if (clusterMinMax[i << 1] <= V && V < clusterMinMax[i << 1 | 1])
        size += clusterSize[i] - 1;
    return size;
  }

  private long possibleSizeGEValue(final long V) {
    long size = 0;
    for (int i = 0; i < clusterNum; i++)
      if (V <= clusterMinMax[i << 1])
        size += clusterSize[i];
      else if (clusterMinMax[i << 1] < V && V <= clusterMinMax[i << 1 | 1])
        size += clusterSize[i] - 1;
    return size;
  }

  public long minValueWithRank(long K) {
    long L = Long.MIN_VALUE, R = Long.MAX_VALUE, mid;
    while (L < R) {
      mid = L + ((R - L) >>> 1);
      if (possibleSizeLEValue(mid) < K) L = mid + 1;
      else R = mid;
    }
//    System.out.println("[minValueWithRank] K:"+K+" val:"+L);
    return L;
  }

  public long maxValueWithRank(long K) {
    K = totN - K + 1;
    long L = Long.MIN_VALUE, R = Long.MAX_VALUE, mid;
    while (L < R) {
      mid = L + ((R - L) >>> 1);
      if (mid == L) mid++;
//      System.out.println("\t\t\t"+L+"  "+mid+"  "+R);
      if (possibleSizeGEValue(mid) < K) R = mid - 1;
      else L = mid;
    }
//    System.out.println("[maxValueWithRank] K:"+K+" val:"+L);
    return L;
  }

  public List<Long> findResultRange(final long K1, final long K2) {
    if(!sorted)
      sortCluster();
    List<Long> result = new ArrayList<>(4);

//    System.out.print("\t\t\t");
//    for(int i=0;i<clusterNum;i++)System.out.print(clusterMinMax[i<<1]+"~"+clusterMinMax[i<<1|1]+":"+clusterSize[i]+"\t");
//    System.out.println();
//    System.out.print("\t\t\t");
//    for(int i=0;i<clusterNum;i++)System.out.print((clusterMinMax[i<<1])+"~"+(clusterMinMax[i<<1|1])+":"+clusterSize[i]+"\t");
//    System.out.println();
    result.add(minValueWithRank(K1));
    result.add(maxValueWithRank(K1));
    result.add(minValueWithRank(K2));
    result.add(maxValueWithRank(K2));
    return result;
  }
  public long getApproxRank(long v){
//    if(!sorted)sortCluster();
//    if(v<calcMid(0))return 0;
//    if(v>calcMid(clusterNum-1))return totN;
//    int l=0,r=clusterNum-1;
//    while(l<r){
//      int mid=(l+r)>>1;
//      if(v<=calcMid(mid))r=mid;
//      else l=mid+1;
//    }
    long rank = 0;
    for(int i=0;i<clusterNum;i++)if(clusterMinMax[i<<1|1]<v)rank+=clusterSize[i];
    else if(clusterMinMax[i<<1]<v)
      rank+=clusterSize[i]*1.0*(v-clusterMinMax[i<<1])/(clusterMinMax[i<<1|1]-clusterMinMax[i<<1]);
    return rank;
  }


  public void reset() {
    clusterNum = 0;
    totN = 0;
  }

  public void compactBeforeSerialization(){
//    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
    int start_compression = this.maxSeriByte/3/8*4;
    compaction(start_compression);
//    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
    while(clusterNum*3*8>maxSeriByte){
      start_compression=start_compression*4/5;
      compaction(start_compression);
//      System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
    }
//    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
//    System.out.println("\t\tcompactBeforeSeri OVER\t");
  }
  public void show(){
    System.out.print("\t\t[DEBUG TDigest]\t"+clusterNum+" items\t");
    for(int i=0;i<clusterNum;i++)System.out.print((clusterMinMax[i<<1])+"~"+clusterMinMax[i<<1|1]+":"+clusterSize[i]+"\t");
    System.out.println();
  }


}