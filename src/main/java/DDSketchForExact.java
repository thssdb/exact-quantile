import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2LongOpenHashMap;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;

public class DDSketchForExact implements Serializable {
  private double alpha;
  private double gamma;
  private double multiplier;
  private int bucket_num_limit;
  private int threshold_for_compression;

  private Int2LongOpenHashMap positive_buckets;
  private Int2DoubleOpenHashMap positive_buckets_content;
  private Int2LongOpenHashMap negative_buckets;
  private Int2DoubleOpenHashMap negative_buckets_content;
  private double positive_collapse_bound,negative_collapse_bound;
  private long zero_count;

  private static final double MIN_POSITIVE_VALUE = 1e-30;
  private static final double COEFFICIENT = 1.5;
  boolean valid_buckets=false;
  Bucket[] buckets;

  // B bucket:  at most threshold_for_compression=1.5*B buckets,    OpenHashMap loadFactor0.75
  // a bucket:  count&lastV(Int2Long,Int2Double)  2*(4+8)=24Byte;
  public DDSketchForExact(double alpha, int bucket_num_limit) {
//        System.out.println(alpha);
    this.alpha = alpha;
    this.bucket_num_limit = Math.max(bucket_num_limit, 2);
    this.threshold_for_compression = (int) (bucket_num_limit * COEFFICIENT);
//    System.out.println("\t\t\t\tcompression:"+threshold_for_compression+"\t\tlimit="+bucket_num_limit);

    this.gamma = 2 * alpha / (1 - alpha) + 1;
    this.multiplier = Math.log(Math.E) / (Math.log1p(gamma - 1));
    this.positive_buckets = new Int2LongOpenHashMap((int) (bucket_num_limit * 0.5));
    this.positive_buckets_content = new Int2DoubleOpenHashMap((int) (bucket_num_limit * 0.5));
    this.negative_buckets = new Int2LongOpenHashMap((int) (bucket_num_limit * 0.5));
    this.negative_buckets_content = new Int2DoubleOpenHashMap((int) (bucket_num_limit * 0.5));
    this.zero_count = 0;
    this.positive_collapse_bound = -Double.MAX_VALUE;
    this.negative_collapse_bound=Double.MAX_VALUE;
  }

  public void insert(double v) {
    valid_buckets=false;
    if (v > MIN_POSITIVE_VALUE) {
      if (v < positive_collapse_bound) {
        v = positive_collapse_bound;
      }
      int i = (int) Math.ceil(Math.log(v) * multiplier);
      positive_buckets.put(i, positive_buckets.getOrDefault(i, 0L) + 1);
      double lastV=positive_buckets_content.getOrDefault(i,v);
      positive_buckets_content.put(i,lastV==v?v:0);
    } else if (v < -MIN_POSITIVE_VALUE) {
      if (v > negative_collapse_bound) {
        v = negative_collapse_bound;
      }
      int i = (int) Math.ceil(Math.log(-v) * multiplier);
      negative_buckets.put(i, negative_buckets.getOrDefault(i, 0L) + 1);
      double lastV=negative_buckets_content.getOrDefault(i,v);
      negative_buckets_content.put(i,lastV==v?v:0);
//      if(v<-0.007032827385454164)System.out.println("^|^\t\tbucketIndex:\t"+i);
    } else {
      zero_count++;
    }
    if(positive_buckets.size() + negative_buckets.size() > threshold_for_compression)
      collapse(bucket_num_limit);
  }

  private void collapse(int limit) {
    int all_exceed = positive_buckets.size() + negative_buckets.size() - limit;
    int nega_exceed = (int)(1.0*negative_buckets.size()/(positive_buckets.size() + negative_buckets.size())*all_exceed);
    int posi_exceed = all_exceed-nega_exceed;
    if(nega_exceed>0&&nega_exceed==negative_buckets.size()){
      nega_exceed--;
      posi_exceed++;
    }
    if(posi_exceed>0&&posi_exceed==positive_buckets.size()){
      posi_exceed--;
      nega_exceed++;
    }

    int[] indices;
    if(nega_exceed>0) {
      indices = negative_buckets.keySet().toIntArray();
      for(int i=0;i<indices.length;i++)indices[i]=-indices[i];
      Arrays.sort(indices);
      for(int i=0;i<indices.length;i++)indices[i]=-indices[i];
      long count = 0;
      for (int i = Math.max(0, indices.length - nega_exceed); i < indices.length; ++i) {
//        if(-Math.pow(gamma, (indices[i]))<-0.007032827385454164)System.out.println("||||||||||!!!!!!!!!!!!!!!!!!!!FAIL COLLAPSE\t\tbucket:"+indices[i]);
        count += negative_buckets.remove(indices[i]);
        negative_buckets_content.remove(indices[i]);
      }
      if (count > 0) {
        int i = indices.length - nega_exceed - 1;
        if (i >= 0) {
          negative_buckets.put(indices[i], negative_buckets.get(indices[i]) + count);
          negative_buckets_content.put(indices[i], 0);
          negative_collapse_bound = -Math.pow(gamma, indices[i]);
//          if(negative_collapse_bound<-0.007032827385454164)System.out.println("||||||||||!!!!!!!!!!!!!!!!!!!!FAIL COLLAPSE  BOUNDDDDDDDDDD:\t"+negative_collapse_bound);
        } else {
          zero_count += count;
          negative_collapse_bound = 0;
//          System.out.println("\t\tWARN!!! collapse all negative.");
        }
      }
    }

    if (posi_exceed > 0) {
      indices = positive_buckets.keySet().toIntArray();
      Arrays.sort(indices);
      long count = 0;
      for (int i = posi_exceed - 1; i >= 0; --i) {
        count += positive_buckets.remove(indices[i]);
        positive_buckets_content.remove(indices[i]);
      }
      positive_buckets.put(indices[posi_exceed], positive_buckets.get(indices[posi_exceed]) + count);
//      if(indices[posi_exceed]>=467)System.out.println("\t[collapse clear].\t"+indices[posi_exceed]+"\t\tR:"+);
      positive_buckets_content.put(indices[posi_exceed], 0);
      positive_collapse_bound = Math.pow(gamma, indices[posi_exceed]);
    }
  }
  static final int DIVIDE_DELTA = 1000000000,DIVIDE_HALF=DIVIDE_DELTA/2;
  private double getL(int index){
    return index>DIVIDE_HALF?Math.pow(gamma, index-DIVIDE_DELTA-1):(index==DIVIDE_HALF?0:-Math.pow(gamma, index));
  }
  private double getR(int index){
    return index>DIVIDE_HALF?Math.pow(gamma, index-DIVIDE_DELTA):(index==DIVIDE_HALF?0:-Math.pow(gamma, index-1));
  }
  private long getCount(int index){
//    System.out.println("\t\t\t\t\t-index="+(-index));
//    System.out.println("\t\t\t\t\t\t\texist"+negative_buckets.containsKey(-index));
    return index>DIVIDE_HALF?positive_buckets.get(index-DIVIDE_DELTA):(index==DIVIDE_HALF?zero_count:negative_buckets.get(index));
  }
  private double getLastV(int index){
//    System.out.println("\t\t\t\t\t-index="+(-index));
    return index>DIVIDE_HALF?positive_buckets_content.get(index-DIVIDE_DELTA):(index==DIVIDE_HALF?0:negative_buckets_content.get(index));
  }




  private void union_buckets() {
    buckets = new Bucket[sketch_size()];
    int i = 0;
    for (Map.Entry<Integer, Long> e : positive_buckets.entrySet()) {
      buckets[i++] = new Bucket(e.getKey()+DIVIDE_DELTA);
    }
    for (Map.Entry<Integer, Long> e : negative_buckets.entrySet()) {
      buckets[i++] = new Bucket(e.getKey());
    }
    if (zero_count > 0) {
      buckets[i] = new Bucket(DIVIDE_HALF);
    }
    Arrays.sort(buckets, Comparator.comparingDouble(o -> (getL(o.bucketIndex))));
    long sum=0;
    for(i=0;i<sketch_size();i++){
      sum+=getCount(buckets[i].bucketIndex);
      buckets[i].prefixSum = sum;
    }
    valid_buckets=true;
  }

  public long total_count() {
    return positive_buckets.values().longStream().sum() + negative_buckets.values().longStream().sum() + zero_count;
  }

  private int find_p_index(Bucket[] buckets, long total_count, double q) {
    double rank = q * (total_count - 1);
    int tmp1 = Integer.highestOneBit(buckets.length);
    int p=-1;
    while(tmp1>0){
      if(p+tmp1<buckets.length&& buckets[p+tmp1].prefixSum<=rank)
        p+=tmp1;
      tmp1/=2;
    }
    return p+1;
  }
  private int find_p_index(Bucket[] buckets, long query_rank) {
    int tmp1 = Integer.highestOneBit(buckets.length);
    int p=-1;
    while(tmp1>0){
      if(p+tmp1<buckets.length&& buckets[p+tmp1].prefixSum<query_rank)
        p+=tmp1;
      tmp1/=2;
    }
    return p+1;
  }
  private int find_p_index_LEQ(Bucket[] buckets, long query_rank) {
    int tmp1 = Integer.highestOneBit(buckets.length);
    int p=-1;
    while(tmp1>0){
      if(p+tmp1<buckets.length&& buckets[p+tmp1].prefixSum<=query_rank)
        p+=tmp1;
      tmp1/=2;
    }
    return p==-1?0:p;
  }
  private int find_p_index_GEQ(Bucket[] buckets, long query_rank) {
    int tmp1 = Integer.highestOneBit(buckets.length);
    int p=-1;
    while(tmp1>0){
      if(p+tmp1<buckets.length&& buckets[p+tmp1].prefixSum<query_rank)
        p+=tmp1;
      tmp1/=2;
    }
    return p+1;
  }

  public double getQuantile(double q) {
    if(!valid_buckets)union_buckets();
    long total_count = total_count();
    Bucket p = buckets[find_p_index(buckets, total_count, q)];
    if (getL(p.bucketIndex) < 0) {
      return 2 * getL(p.bucketIndex) / (1 + gamma);
    } else {
      return 2 * getR(p.bucketIndex) / (1 + gamma);
    }
  }


  public double[] findResultRange(long K1, long K2) {
    if(!valid_buckets)union_buckets();
    DoubleArrayList result = new DoubleArrayList(2);
    int p1=find_p_index/*_LEQ*/(buckets,K1),p2=find_p_index/*_GEQ*/(buckets,K2);
//    System.out.println("\t\tfinding range\t\tK1,2:"+K1+","+K2+"\t\tp1,2:"+p1+","+p2+"\t\tprefixSum1,2:"+buckets[p1].prefixSum+","+buckets[p2].prefixSum);
//    System.out.println("\t\t\t"+"num_nega,posi_bucket:"+negative_buckets.size()+","+positive_buckets.size());
//    System.out.println("\t\t\t"+"nega,posi_collapse_bound:"+negative_collapse_bound+","+positive_collapse_bound);
//    System.out.println("\t\t\tbucket1\tindex,L,R:"+buckets[p1].bucketIndex+","+getL(buckets[p1].bucketIndex)+","+getR(buckets[p1].bucketIndex));
//    System.out.println("\t\t\tbucket2\tindex,L,R:"+buckets[p2].bucketIndex+","+getL(buckets[p2].bucketIndex)+","+getR(buckets[p2].bucketIndex));

//    int tmpCount=0;
//    for(int i=0;i<p1;i++) {
//      tmpCount+=getCount(buckets[i].bucketIndex);
//      System.out.println("\t\t\t\tMMP small buckets.\tL:" + getL(buckets[i].bucketIndex)+"\tCount:\t"+getCount(buckets[i].bucketIndex)+"\t\tindex:"+buckets[i].bucketIndex);
//    }System.out.println("\t\t\t\tSmallerCount:"+tmpCount);

    double valL=getL(buckets[p1].bucketIndex),valR=getR(buckets[p2].bucketIndex);
    if(valL>0&&valL<=positive_collapse_bound)valL=MIN_POSITIVE_VALUE;
    if(valR<0&&valR>=negative_collapse_bound)valR=-MIN_POSITIVE_VALUE;
//    if(valL==0)valL=p1==0?-Double.MAX_VALUE:(getR(buckets[p1-1].bucketIndex)+Double.MIN_NORMAL);
//    if(valR==0)valL=(p2==buckets.length-1)?Double.MAX_VALUE:(getL(buckets[p2+1].bucketIndex)-Double.MIN_NORMAL);
    if(valL==0&&valR==0){
      result.add(0);
      result.add(0);
      result.add(-233);
    }else {
      if (valL == 0) valL = -MIN_POSITIVE_VALUE;
      if (valR == 0) valR = MIN_POSITIVE_VALUE;
      if (getLastV(buckets[p1].bucketIndex) != 0 && getLastV(buckets[p2].bucketIndex) != 0) {
//        System.out.println("\t\tWIN!!!");
        result.add(getLastV(buckets[p1].bucketIndex));
        result.add(getLastV(buckets[p2].bucketIndex));
        result.add(-233);
      }else {
        result.add(valL);
        result.add(valR);
      }
    }
    return result.toDoubleArray();
  }


//  public double[] getFilterL(long CountOfValL,long CountOfValR,double valL,double valR,long K){
//    if(K<=CountOfValL)return new double[]{valL,-233.0};
//    if(K>CountOfValL+total_count())return new double[]{valR,-233.0};
//    K-=CountOfValL;
//    int p1=find_p_index(buckets,K);
//    valL=getL(buckets[p1].bucketIndex);
//    if(valL>0&&valL<=positive_collapse_bound)valL=MIN_POSITIVE_VALUE;
//    if(valL==0)return new double[]{0,-233};
//    if(getLastV(buckets[p1].bucketIndex) != 0)return new double[]{getLastV(buckets[p1].bucketIndex),-233};
//    return new double[]{valL};
//  }
//
//  public double[] getFilter(long CountOfValL,long CountOfValR,double lastValL,double lastValR,long K1,long K2) {
//    if(!valid_buckets)union_buckets();
//    DoubleArrayList result = new DoubleArrayList(2);
//    int p1=find_p_index/*_LEQ*/(buckets,K1),p2=find_p_index/*_GEQ*/(buckets,K2);
//
//    double valL=getL(buckets[p1].bucketIndex),valR=getR(buckets[p2].bucketIndex);
//    if(valL>0&&valL<=positive_collapse_bound)valL=MIN_POSITIVE_VALUE;
//    if(valR<0&&valR>=negative_collapse_bound)valR=-MIN_POSITIVE_VALUE;
//    if(valL==0&&valR==0){
//      result.add(0);
//      result.add(0);
//      result.add(-233);
//    }else {
//      if (valL == 0) valL = -MIN_POSITIVE_VALUE;
//      if (valR == 0) valR = MIN_POSITIVE_VALUE;
//      if (getLastV(buckets[p1].bucketIndex) != 0 && getLastV(buckets[p2].bucketIndex) != 0) {
////        System.out.println("\t\tWIN!!!");
//        result.add(getLastV(buckets[p1].bucketIndex));
//        result.add(getLastV(buckets[p2].bucketIndex));
//        result.add(-233);
//      }else {
//        result.add(valL);
//        result.add(valR);
//      }
//    }
//    return result.toDoubleArray();
//  }

  public int sketch_size() {
    return positive_buckets.size() + negative_buckets.size() + (zero_count == 0 ? 0 : 1);
  }

  private static class Bucket {
    public int bucketIndex;
    public long prefixSum;

    Bucket(int bucketIndex) {
      this.bucketIndex=bucketIndex;
      this.prefixSum=0;
    }
  }
}
