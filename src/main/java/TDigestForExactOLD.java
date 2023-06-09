//// inspired by t-Digest by Ted Dunning. See https://github.com/tdunning/t-digest
//// This is a simple implementation with radix sort and K0.
//// Clusters are NOT strictly in order.
//
//import it.unimi.dsi.fastutil.ints.IntArrayList;
//import org.apache.iotdb.tsfile.utils.ReadWriteIOUtils;
//import org.eclipse.collections.api.tuple.primitive.DoubleLongPair;
//import org.eclipse.collections.impl.tuple.primitive.PrimitiveTuples;
//
//import java.io.IOException;
//import java.io.InputStream;
//import java.io.OutputStream;
//import java.nio.ByteBuffer;
//import java.text.DecimalFormat;
//import java.util.List;
//
//public class TDigestForExactOLD {
//  public final int compression;
////  public DoubleLongPair[] cluster;
//  public double[] clusterAvgCountMinMax;
//  IntArrayList index;
//  public int clusterNum, clusterNumMemLimit, clusterNumSeriLimit;
//  public int maxSeriByte, maxMemByte;
//  public long totN;
//  boolean sorted = true;
//
////  public TDigestForExact(int maxMemByte, int maxSeriByte) {
////    this.maxMemByte = maxMemByte;
////    this.maxSeriByte = maxSeriByte;
////    this.clusterNumMemLimit = maxMemByte / 5 / 8;
////    this.clusterNumSeriLimit = maxSeriByte / (8+2);
////    this.compression = this.clusterNumMemLimit / 8; //  cluster:buffer = 1:7
////    clusterNum = 0;
////    totN = 0;
////    cluster = new DoubleLongPair[clusterNumMemLimit];
////    cluster = new DoubleLongPair[clusterNumMemLimit];
////  }
//  public TDigestForExactOLD(int maxMemByte) {
//    this.maxMemByte = maxMemByte;
//    this.maxSeriByte = maxMemByte;
//    this.clusterNumMemLimit = maxMemByte / 8 / 4;
//    this.clusterNumSeriLimit = maxSeriByte / 2 / 8;
//    this.compression = this.clusterNumMemLimit / 8; //  cluster:buffer = 1:7
//    clusterNum = 0;
//    totN = 0;
//    clusterAvgCountMinMax = new double[clusterNumMemLimit<<2];
//    index = new IntArrayList(clusterNumMemLimit);
//  }
//
//  private void addCluster(double avg,long count,double min,double max) {
//    if (clusterNum == clusterNumMemLimit) compaction(this.compression);
//    clusterAvgCountMinMax[clusterNum<<2] = avg;
//    clusterAvgCountMinMax[clusterNum<<2|1] = count;
//    clusterAvgCountMinMax[clusterNum<<2|2] = min;
//    clusterAvgCountMinMax[clusterNum<<2|3] = max;
//    index.add(clusterNum++);
//  }
//  private void addCluster(double v) {
//    if (clusterNum == clusterNumMemLimit) compaction(this.compression);
//    clusterAvgCountMinMax[clusterNum<<2] = v;
//    clusterAvgCountMinMax[clusterNum<<2|1] = 1;
//    clusterAvgCountMinMax[clusterNum<<2|2] = v;
//    clusterAvgCountMinMax[clusterNum<<2|3] = v;
//    index.add(clusterNum++);
//  }
//
//  DoubleLongPair mergeTwoCluster(DoubleLongPair a,DoubleLongPair b){
//    long sizeA = a.getTwo(),sizeB = b.getTwo();
//    return PrimitiveTuples.pair(a.getOne()*(1.0*sizeA/(sizeA+sizeB))+b.getOne()*(1.0*sizeB/(sizeA+sizeB)),a.getTwo()+b.getTwo());
//  }
//  double mergeAvg(double avg1,long count1,double avg2,long count2){
//    return avg1*(1.0*count1/(count1+count2))+avg2*(1.0*count2/(count1+count2));
//  }
//
//  public void sortCluster() {
//    if(sorted)return;
//    index.sort((x,y)->(Double.compare(clusterAvgCountMinMax[x<<2],clusterAvgCountMinMax[y<<2])));
//    sorted = true;
//  }
//  // TODO
//  private void compaction(int compression) {
//
//    sortCluster();
//    int tmpClusterNum = clusterNum;
//    clusterNum = 0;
//
//    long expectedClusterSize = (totN + compression - 1) / compression;
//    DoubleLongPair cnt = PrimitiveTuples.pair(0d, 0L);
////    System.out.println("\t\t expectedClusterSize:"+expectedClusterSize + " totN:"+totN+" oldClusterNum:"+tmpClusterNum);
//    for (int i = 0; i < tmpClusterNum; i++) {
//      if(compression!=this.compression&&tmpClusterNum-i+clusterNum+(cnt.getTwo()>0?1:0)<=maxSeriByte/2/8){
////        System.out.println("\t\t\t\tpartial compaction!");
//        if(cnt.getTwo()>0)
//          addCluster(cnt);
//        for(;i<tmpClusterNum;i++)
//          addCluster(cluster[i]);
////        System.out.println("\t\t\t\tpartial compaction!"+clusterNum);
//        return;
//      }
//      if (cnt.getTwo() + cluster[i].getTwo() <= expectedClusterSize) {
//        cnt = mergeTwoCluster(cnt,cluster[i]);
//      } else {
//        if (cnt.getTwo() > 0)
//          addCluster(cnt);
//        cnt = cluster[i];
//      }
//    }
//    if(cnt.getTwo()>0)
//      addCluster(cnt);
////    System.out.println("\t\t after compaction:"+tmpClusterNum+"-->"+clusterNum);
////    show();
//  }
//
//
////  public void add(final long value,final long count){
//////    System.out.println("\t\t add: "+(value)+"  "+count);
////    totN += count;
//////    bufferValueCount[bufferNum<<1]=value;
//////    bufferValueCount[bufferNum<<1|1]=count;
////    bufferNum++;
////    if(bufferNum == bufferSizeLimit)
////      updateFromBuffer();
////  }
//
////  public void update(long value) {
//////    System.out.println("\t\t add: "+(value));
////    totN++;
//////    value ^= deltaForUnsigned;
////    sorted = false;
////    addCluster(1, value, value);
////  }
//  public void update(double value) {
////    System.out.println("\t\t add: "+(value));
//    totN++;
//    sorted = false;
//    addCluster(value);
//  }
//
//  public void merge(TDigestForExactOLD another) {
//    sorted = false;
////    System.out.println("\t\t add: "+(value));
//    totN += another.totN;
//    for (int i = 0; i < another.clusterNum; i++)
//      addCluster(another.cluster[i]);
//  }
//  public void merge(List<TDigestForExactOLD> anotherList) {
//    for(TDigestForExactOLD another:anotherList)
//      merge(another);
//  }
//  public double quantile(double q){
//    sortCluster();
//    double preN=0;
//    for(int i=0;i<clusterNum;i++) {
//
//      if(preN+0.5*cluster[i].getTwo()>=q*totN) {
//        if (i == 0) return cluster[i].getOne();
//        DoubleLongPair c1 = cluster[i-1],c2=cluster[i];
//        double wLeft = q * totN - preN + 0.5 * c1.getTwo();
//        double wRight = preN - q * totN + 0.5 * c2.getTwo();
//        return (c1.getOne()*wRight+c2.getOne()*wLeft)/(wLeft+wRight);
//      }
//      preN += cluster[i].getTwo();
//    }
//    return cluster[clusterNum-1].getOne();
//  }
//
//
//
//
//  public void reset() {
//    clusterNum = 0;
//    totN = 0;
//  }
//
//  public void compactBeforeSerialization(){
////    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
//    int start_compression = this.maxSeriByte/(8+2)*4;
//    compaction(start_compression);
////    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
//    while(clusterNum*2*8>maxSeriByte){
//      start_compression=start_compression*4/5;
////      System.out.println("\t\t?compactBeforeSeri\t\t\t\ttry compression"+start_compression);
//      compaction(start_compression);
////      System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
//    }
////    System.out.println("\t\t?compactBeforeSeri\t"+clusterNum);
////    System.out.println("\t\tcompactBeforeSeri OVER\t");
//  }
//  public void show(){
//    System.out.print("\t\t[DEBUG TDigest]\t"+clusterNum+" items\t");
//    DecimalFormat df = new DecimalFormat("0.0E000");
//    for(int i=0;i<clusterNum;i++)System.out.print("("+ df.format(cluster[i].getOne())+","+cluster[i].getTwo()+")"+"\t");
//    System.out.println();
//  }
//
//  //
//
//  public int serialize(OutputStream outputStream) throws IOException {// 15+1*?+8*?
//    compactBeforeSerialization();// if N==maxN
//    int byteLen = 0;
//    byteLen+=ReadWriteIOUtils.write(totN, outputStream);
//    byteLen+=ReadWriteIOUtils.write((short)clusterNum, outputStream);
//    for(int i=0;i<clusterNum;i++){
//      byteLen+=ReadWriteIOUtils.write(cluster[i].getOne(),outputStream);
//      byteLen+=ReadWriteIOUtils.write((short)cluster[i].getTwo(),outputStream);
//    }
//    return byteLen;
//  }
//
//  public TDigestForExactOLD(InputStream inputStream, int maxMemoryByte, int maxSerializeByte) throws IOException {
//    this(maxMemoryByte);
//    this.totN = ReadWriteIOUtils.readLong(inputStream);
//    int clusterNum = ReadWriteIOUtils.readShort(inputStream);
//    for(int i=0;i<clusterNum;i++){
//      double a = ReadWriteIOUtils.readDouble(inputStream);
//      long b = ReadWriteIOUtils.readShort(inputStream);
//      addCluster(PrimitiveTuples.pair(a,b));
//    }
//    this.sorted = false;
//  }
//
//  public TDigestForExactOLD(ByteBuffer byteBuffer, int maxMemoryByte, int maxSerializeByte) {
//    this(maxMemoryByte);
//    this.totN = ReadWriteIOUtils.readLong(byteBuffer);
//    int clusterNum = ReadWriteIOUtils.readShort(byteBuffer);
//    for(int i=0;i<clusterNum;i++){
//      double a = ReadWriteIOUtils.readDouble(byteBuffer);
//      long b = ReadWriteIOUtils.readShort(byteBuffer);
//      addCluster(PrimitiveTuples.pair(a,b));
//    }
//    this.sorted = false;
//  }
//
//  public TDigestForExactOLD(InputStream inputStream) throws IOException {
//    this.totN = ReadWriteIOUtils.readLong(inputStream);
//    int clusterNum = ReadWriteIOUtils.readShort(inputStream);
//    int maxMemByte = clusterNum*16+16;
//    this.maxMemByte = maxMemByte;
//    this.maxSeriByte = maxMemByte;
//    this.clusterNumMemLimit = maxMemByte / 2 / 8;
//    this.clusterNumSeriLimit = maxSeriByte / (2 + 8);
////    this.compression = this.clusterNumMemLimit / 6; //  cluster:buffer = 1:5
//    this.compression = this.clusterNumMemLimit / 8; //  cluster:buffer = 1:7
//    this.clusterNum = 0;
//    cluster = new DoubleLongPair[clusterNumMemLimit];
//    for(int i=0;i<clusterNum;i++){
//      double a = ReadWriteIOUtils.readDouble(inputStream);
//      long b = ReadWriteIOUtils.readShort(inputStream);
//      addCluster(PrimitiveTuples.pair(a,b));
//    }
//    this.sorted = false;
//  }
//
//  public TDigestForExactOLD(ByteBuffer byteBuffer) {
//
//    this.totN = ReadWriteIOUtils.readLong(byteBuffer);
//    int clusterNum = ReadWriteIOUtils.readShort(byteBuffer);
//    int maxMemByte = clusterNum * 16 + 16;
//    this.maxMemByte = maxMemByte;
//    this.maxSeriByte = maxMemByte;
//    this.clusterNumMemLimit = maxMemByte / 2 / 8;
//    this.clusterNumSeriLimit = maxSeriByte / (2 + 8);
////    this.compression = this.clusterNumMemLimit / 6; //  cluster:buffer = 1:5
//    this.compression = this.clusterNumMemLimit / 8; //  cluster:buffer = 1:7
//    this.clusterNum = 0;
//    cluster = new DoubleLongPair[clusterNumMemLimit];
//    for (int i = 0; i < clusterNum; i++) {
//      double a = ReadWriteIOUtils.readDouble(byteBuffer);
//      long b = ReadWriteIOUtils.readShort(byteBuffer);
//      addCluster(PrimitiveTuples.pair(a, b));
//    }
//    this.sorted = false;
//  }
//}