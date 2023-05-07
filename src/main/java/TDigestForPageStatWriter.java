import org.apache.iotdb.tsfile.utils.ReadWriteIOUtils;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;

import java.io.IOException;
import java.io.OutputStream;

public class TDigestForPageStatWriter extends TDigestForPageStat{
  int N,maxN;
  short stat_K,rnd; // rank of first num
  LongArrayList num;


  public TDigestForPageStatWriter(int pageMaxNum,short K) { // K is a power of 2
    maxN = pageMaxNum;
    num = new LongArrayList(maxN);
    stat_K=K;
  }

  public int getN(){return N;}

  public void update(long x){ // signed long
    N++;
    num.add(x);
  }

  public int serialize(OutputStream outputStream) throws IOException {
    rnd = (short) (getNextRand(stat_K-1)+1); // rnd=1..stat_K
    num.sortThis();
    int byteLen = 0;
    byteLen+= ReadWriteIOUtils.write(N, outputStream);
    byteLen+= ReadWriteIOUtils.write(stat_K, outputStream);
    byteLen+= ReadWriteIOUtils.write(rnd, outputStream);

    long L,R,cntV,interval_v;
    int rankL=0,rankR = Math.min(N,rnd)-1;
    int cntB=0;
    int debug_n = 0;
    byteLen+= ReadWriteIOUtils.write(num.get(0), outputStream);
//    System.out.println("\t{"+num.get(0)+"}");
//    debug_n++;
    if(rankR==0)
      rankR = Math.min(stat_K,N-1);
    while(rankL<rankR){
      L = num.get(rankL);
      R = num.get(rankR);
      byteLen+= ReadWriteIOUtils.write(R, outputStream);
//      System.out.println("\t{"+R+"}");
//      debug_n++;
      interval_v = L==R?1:(((R-L)>>>BUCKET_K)+1);
      for(int i=rankL+1;i<rankR;i++){ // rank in (rankL,rankR)
        cntV = num.get(i);
        while(cntV>=L+interval_v){  // value in [L,L+interval_v)
          byteLen+= ReadWriteIOUtils.write((byte)cntB, outputStream);
//          System.out.println("\t\t["+L+","+(L+interval_v-1)+"]"+":"+cntB);
//          debug_n+=cntB;
          cntB=0;
          L+=interval_v;
        }
        cntB++;
      }
      while(L<=R){
//        System.out.println("\t\t\t!!!!\t\t\t");
        byteLen+= ReadWriteIOUtils.write((byte)cntB, outputStream);
//        System.out.println("\t\t__["+L+","+((L<=R-interval_v)?(L+interval_v-1):R)+"]"+":"+cntB);
//        debug_n+=cntB;
        cntB=0;
        L+=interval_v;
      }
      rankL=rankR;
      rankR = Math.min(rankL+stat_K,N-1);
    }
//    System.out.println("\t\tcheck n: "+debug_n+" == "+N);

    return byteLen;
  }


}