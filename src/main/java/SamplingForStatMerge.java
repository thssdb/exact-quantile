// inspired by t-Digest by Ted Dunning. See https://github.com/tdunning/t-digest
// This is a simple implementation with radix sort and K0.
// Clusters are NOT strictly in order.

import it.unimi.dsi.fastutil.doubles.Double2ReferenceAVLTreeMap;
import it.unimi.dsi.fastutil.doubles.Double2ReferenceMap;
import it.unimi.dsi.fastutil.objects.ReferenceArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.eclipse.collections.api.tuple.primitive.LongDoublePair;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;
import org.eclipse.collections.impl.tuple.primitive.PrimitiveTuples;

import java.util.Comparator;
import java.util.List;
//"Weighted random sampling with a reservoir"
public class SamplingForStatMerge {
  public final int sampleLimit;
  public Double2ReferenceAVLTreeMap<LongDoublePair> sampleTreeMap;
  public int maxSeriByte, maxMemByte;
  public long totN;
  public XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom();
  boolean sorted = false;
  public ReferenceArrayList<LongDoublePair> sampleList;

  public SamplingForStatMerge(int maxMemByte, int maxSeriByte) {
    this.maxMemByte = maxMemByte;
    this.maxSeriByte = maxSeriByte;
    this.sampleLimit = maxSeriByte/(8);
    this.sampleTreeMap = new Double2ReferenceAVLTreeMap<>();
  }
  public SamplingForStatMerge(int maxMemByte) {
    this.maxMemByte = maxMemByte;
    this.maxSeriByte = maxMemByte;
    this.sampleLimit = maxMemByte/(8*3);
    this.sampleTreeMap = new Double2ReferenceAVLTreeMap<>();
  }

  private void putSample(double score, LongDoublePair pair){
    LongDoublePair old = sampleTreeMap.put(score, pair);
    for(int i=1;(old!=null);i++){
//      System.out.println("\t\t??");
      score += (((i&1)==0)?Double.MIN_NORMAL:-Double.MIN_NORMAL)*i*i;
      old = sampleTreeMap.put(score, old);
    }
  }

  private void add(long value,double weight){
    double score = weight==1?random.nextDoubleFast():Math.pow(random.nextDoubleFast(),1.0/weight);
    if(sampleTreeMap.size()<sampleLimit)putSample(score, PrimitiveTuples.pair(value,weight));
    else if(sampleTreeMap.firstDoubleKey()<score){
      sampleTreeMap.remove(sampleTreeMap.firstDoubleKey());
      putSample(score,PrimitiveTuples.pair(value,weight));
    }
  }

  public void update(long value) {
    totN++;
    add(value,1);
  }

  public void merge(SamplingForStatMerge another) {
    totN += another.totN;
    if(another.sampleTreeMap!=null) {
      for (LongDoublePair pair : another.sampleTreeMap.values())
        add(pair.getOne(), pair.getTwo());
    }else{
      for (LongDoublePair pair : another.sampleList)
        add(pair.getOne(), pair.getTwo());
    }
  }
  public void merge(List<SamplingForStatMerge> anotherList) {
    for(SamplingForStatMerge another:anotherList)
      merge(another);
  }
  public void sortSample(){
    if(sorted)return;
    if(sampleTreeMap!=null)
      sampleList = new ReferenceArrayList<>(sampleTreeMap.values());
    sampleList.sort(Comparator.comparing(LongDoublePair::getOne));
    sorted=true;
    sampleTreeMap=null;
  }

  public long quantile(double q){
    sortSample();
    if(q<=0)return sampleList.get(0).getOne();
    if(q>=1)return sampleList.top().getOne();
    return sampleList.get((int)(q*sampleList.size())).getOne();
  }




  public void reset() {
    totN = 0;
    sorted = false;
    sampleList = null;
    this.sampleTreeMap = new Double2ReferenceAVLTreeMap<>();
  }

  public void compactBeforeSerialization(){
    // no-op
    double weight = 1.0*totN/sampleTreeMap.size();
    sampleList = new ReferenceArrayList<>(sampleTreeMap.size());
//    sampleTreeMap.forEach((x,y)->(sampleList.add(PrimitiveTuples.pair(y.getOne(),weight))));
    for(Double2ReferenceMap.Entry<LongDoublePair> entry:sampleTreeMap.double2ReferenceEntrySet())
      sampleList.add(PrimitiveTuples.pair(entry.getValue().getOne(),weight));
    sampleTreeMap = null;
  }
  public void show(){
    sortSample();System.out.print("\t\t["+sampleList.size()+" samples for N="+totN+"]");
    for(int i=0;i<sampleList.size();i++)System.out.print("\t("+sampleList.get(i).getOne()+","+sampleList.get(i).getTwo()+")");
    System.out.println();
  }
  public void showNum(){
    sortSample();System.out.println("\t\t["+sampleList.size()+" samples for N="+totN+"]");
  }
  public int getSampleSize(){
    return sampleTreeMap==null?sampleList.size():sampleTreeMap.size();
  }
  public double getAvgErr(){
    return 1/3.0*Math.pow(getSampleSize(),-0.5);
  }
  public LongArrayList getLowerBound(double q){
    LongArrayList lowerBound = new LongArrayList();
    for(int k=4;k<=6;k+=2)
      lowerBound.add(quantile(q-k*getAvgErr()));
    return lowerBound;
  }
  public LongArrayList getUpperBound(double q){
    LongArrayList upperBound = new LongArrayList();
    for(int k=4;k<=6;k+=2)
      upperBound.add(quantile(q+k*getAvgErr()));
    return upperBound;
  }



}