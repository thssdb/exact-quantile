// inspired by t-Digest by Ted Dunning. See https://github.com/tdunning/t-digest
// This is a simple implementation with radix sort and K0.
// Clusters are NOT strictly in order.

import it.unimi.dsi.fastutil.objects.ReferenceArrayList;
import org.eclipse.collections.api.tuple.primitive.DoubleLongPair;
import org.eclipse.collections.impl.tuple.primitive.PrimitiveTuples;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
// 没用。
public class TDigestListForStatMerge {
  public final long deltaForUnsigned;
  public final int compression;public ReferenceArrayList<DoubleLongPair> cluster;
  public int clusterNum, clusterNumMemLimit, clusterNumSeriLimit;
  public int maxSeriByte, maxMemByte;
  public long totN;
  boolean sorted = true;

  public TDigestListForStatMerge(int maxMemByte, int maxSeriByte) {
    this.maxMemByte = maxMemByte;
    this.maxSeriByte = maxSeriByte;
    this.clusterNumMemLimit = maxMemByte / 2 / 8;
    clusterNumSeriLimit = maxSeriByte / 2 / 8;
    this.compression = this.clusterNumMemLimit / 6; //  cluster:buffer = 1:5
    clusterNum = 0;
    totN = 0;
    cluster = new ReferenceArrayList<>(clusterNumMemLimit);
    cluster.size(clusterNumMemLimit);
    deltaForUnsigned = 1L << (64 - 1);
  }

  private void addCluster(DoubleLongPair pair) {
    cluster.set(clusterNum,pair);
    clusterNum++;
    if (clusterNum == clusterNumMemLimit) compaction(this.compression);
  }
  private void addCluster(double v){addCluster(PrimitiveTuples.pair(v, 1L));}
  DoubleLongPair mergeTwoCluster(DoubleLongPair a,DoubleLongPair b){
    long sizeA = a.getTwo(),sizeB = b.getTwo();
    return PrimitiveTuples.pair(a.getOne()*(1.0*sizeA/(sizeA+sizeB))+b.getOne()*(1.0*sizeB/(sizeA+sizeB)),a.getTwo()+b.getTwo());
  }

  public void sortCluster() {
    if(sorted)return;
    cluster.size(clusterNum);
    cluster.sort(Comparator.comparingDouble(DoubleLongPair::getOne));
    cluster.size(clusterNumMemLimit);
    sorted = true;
  }

  private void compaction(int compression) {

    sortCluster();
    int tmpClusterNum = clusterNum;
    clusterNum = 0;

    long expectedClusterSize = (totN + compression - 1) / compression;
    DoubleLongPair cnt = PrimitiveTuples.pair(0d, 0L);
//    System.out.println("\t\t expectedClusterSize:"+expectedClusterSize + " totN:"+totN+" oldClusterNum:"+tmpClusterNum);
    for (int i = 0; i < tmpClusterNum; i++) {
      if(compression!=this.compression&&tmpClusterNum-i+clusterNum+(cnt.getTwo()>0?1:0)<=maxSeriByte/2/8){
//        System.out.println("\t\t\t\tpartial compaction!");
        if(cnt.getTwo()>0)
          addCluster(cnt);
        for(;i<tmpClusterNum;i++)
          addCluster(cluster.get(i));
//        System.out.println("\t\t\t\tpartial compaction!"+clusterNum);
        return;
      }
      if (cnt.getTwo() + cluster.get(i).getTwo() <= expectedClusterSize) {
        cnt = mergeTwoCluster(cnt,cluster.get(i));
      } else {
        if (cnt.getTwo() > 0)
          addCluster(cnt);
        cnt = cluster.get(i);
      }
    }
    if(cnt.getTwo()>0)
      addCluster(cnt);
//    System.out.println("\t\t after compaction:"+tmpClusterNum+"-->"+clusterNum);
//    show();
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

//  public void update(long value) {
////    System.out.println("\t\t add: "+(value));
//    totN++;
////    value ^= deltaForUnsigned;
//    sorted = false;
//    addCluster(1, value, value);
//  }
  public void update(double value) {
//    System.out.println("\t\t add: "+(value));
    totN++;
    sorted = false;
    addCluster(value);
  }

  public void merge(TDigestListForStatMerge another) {
    sorted = false;
//    System.out.println("\t\t add: "+(value));
    totN += another.totN;
    for (int i = 0; i < another.clusterNum; i++)
      addCluster(another.cluster.get(i));
  }
  public void merge(List<TDigestListForStatMerge> anotherList) {
    for(TDigestListForStatMerge another:anotherList)
      merge(another);
  }
  public double quantile(double q){
    sortCluster();
    double preN=0;
    for(int i=0;i<clusterNum;i++) {

      if(preN+0.5*cluster.get(i).getTwo()>=q*totN) {
        if (i == 0) return cluster.get(i).getOne();
        DoubleLongPair c1 = cluster.get(i-1),c2=cluster.get(i);
        double wLeft = q * totN - preN + 0.5 * c1.getTwo();
        double wRight = preN - q * totN + 0.5 * c2.getTwo();
        return (c1.getOne()*wRight+c2.getOne()*wLeft)/(wLeft+wRight);
      }
      preN += cluster.get(i).getTwo();
    }
    return cluster.get(clusterNum-1).getOne();
  }




  public void reset() {
    clusterNum = 0;
    totN = 0;
  }

  public void compactBeforeSerialization(){
//    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
    int start_compression = this.maxSeriByte/2/8*4;
    compaction(start_compression);
//    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
    while(clusterNum*2*8>maxSeriByte){
      start_compression=start_compression*4/5;
//      System.out.println("\t\t?compactBeforeSeri\t\t\t\ttry compression"+start_compression);
      compaction(start_compression);
//      System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
    }
//    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
//    System.out.println("\t\tcompactBeforeSeri OVER\t");
  }
  public void show(){
    System.out.print("\t\t[DEBUG TDigest]\t"+clusterNum+" items\t");
    DecimalFormat df = new DecimalFormat("0.0E000");
    for(int i=0;i<clusterNum;i++)System.out.print("("+ df.format(cluster.get(i).getOne())+","+cluster.get(i).getTwo()+")"+"\t");
    System.out.println();
  }


}