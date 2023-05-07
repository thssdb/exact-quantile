import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class SamplerFastForMedian {
//  long[] sampleValueCount;
  long[] sample;
  final int sampleSizeLimit;
  int sampleNum;
  long totSize;
  long randSeed;
  long minValue,maxValue;
  Random random = new Random(2333);

  public SamplerFastForMedian(int sampleSizeLimit) {
    sample = new long[sampleSizeLimit];
    this.sampleSizeLimit = sampleSizeLimit;
    reset();
  }
  private long nextLong(){
    return random.nextLong();
//    randSeed ^= (randSeed << 21);
//    randSeed ^= (randSeed >>> 35);
//    randSeed ^= (randSeed << 4);
//    return randSeed;
  }


  public void add(final long value){
    minValue = Math.min(minValue, value);
    maxValue = Math.max(maxValue, value);
    totSize++;
    if(sampleNum<sampleSizeLimit)
      sample[sampleNum++] = value;
    else {
      long p = Math.abs(nextLong())%totSize;
      if(p<sampleSizeLimit)
        sample[(int)p] = value;
    }
  }
  public boolean isExactResult(){
    return sampleNum==totSize;
  }

  public List<Long> findResultRange(long K1, long K2) {
    List<Long> result = new ArrayList<>(2);
    Arrays.sort(sample, 0, sampleNum);
    if(isExactResult()){
      result.add(sample[(int)(K1-1)]);
      result.add(sample[(int)(K2-1)]);
      return result;
    }
    double weight = 1.0d * totSize / sampleNum;
    int posL1 = (int) Math.floor((K1-1) / weight)-1, posR2 = (int) Math.ceil(K2 / weight);
//    int delta = (int)Math.ceil(Math.pow(Math.log(sampleNum),2.0))*8;
    int delta = (int)Math.pow(sampleNum,0.33)+100;
    posL1-=delta;posR2+=delta;
    if(posL1<0)
      result.add(minValue);
    else
      result.add(sample[posL1]);
    if(posR2>=sampleNum)
      result.add(maxValue);
    else
      result.add(sample[posR2]);
    System.out.println("\t\t\t\t sampler find: w:"+weight+"  sampleSize:"+sampleNum+"  posL1:"+posL1+"  posR2:"+posR2+"    K1,K2:"+K1+","+K2);
    System.out.println("\t\t\t\t ");
    return result;
  }
  public void reset() {
    sampleNum = 0;
    totSize = 0;
    randSeed = 233L << 16L;
    minValue = Long.MAX_VALUE;
    maxValue = Long.MIN_VALUE;
  }

}