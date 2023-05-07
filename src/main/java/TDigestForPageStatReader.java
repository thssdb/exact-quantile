import org.apache.iotdb.tsfile.utils.ReadWriteIOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;

public class TDigestForPageStatReader extends TDigestForPageStat{
  short stat_K,rnd; // rank of first num
  int N;
  long[] clusterMinMaxSize;
  int clusterNum;
  int readIndex;

  private void addCluster(long L,long R,long times){
    if(times==0)return;
    int index=(clusterNum++)*3;
    clusterMinMaxSize[index++]=L;
    clusterMinMaxSize[index++]=R;
    clusterMinMaxSize[index]=times;
  }
  private void addSingleNum(long num){
    int index=(clusterNum++)*3;
    clusterMinMaxSize[index++]=num;
    clusterMinMaxSize[index++]=num;
    clusterMinMaxSize[index]=1;
  }

  public boolean hasNextCluster(){return readIndex<clusterNum*3;}
  public long readNext(){return clusterMinMaxSize[readIndex++];}
  public boolean hasPossibleValueIn(long L,long R){
    boolean flag=false;
    for(int i=0;i<clusterNum*3&&!flag;i+=3)
      flag = !(clusterMinMaxSize[i] < L || clusterMinMaxSize[i + 1] > R);
    return flag;
  }


  public TDigestForPageStatReader(InputStream inputStream) throws IOException {
    this.N = ReadWriteIOUtils.readInt(inputStream);
    this.stat_K = ReadWriteIOUtils.readShort(inputStream);
    this.rnd = ReadWriteIOUtils.readShort(inputStream);// rnd=1..stat_K
    int bucketApproxNum = (N/stat_K+3);
    int clusterApproxNum = bucketApproxNum*(1+(1<<BUCKET_K));
    clusterMinMaxSize = new long[clusterApproxNum*3+5];
    clusterNum=readIndex=0;
//    System.out.println("\t\treading TDigest.\tN:"+N+"\tstat_K:"+stat_K+"\trnd:"+rnd);
    long L,R,interval_v;
    int rankL=0,rankR = Math.min(N,rnd)-1;
    int cntB;
    int debug_n = 0;
    L = ReadWriteIOUtils.readLong(inputStream);
    addSingleNum(L);
//    System.out.println("\t{"+L+"}");
//    debug_n++;
    if(rankR==0)
      rankR = Math.min(stat_K,N-1);
    while(rankL<rankR){
      R = ReadWriteIOUtils.readLong(inputStream);
//      System.out.println("\t{"+R+"}");
      debug_n++;
      interval_v = L==R?1:(((R-L)>>>BUCKET_K)+1);
      while(L<=R){
        cntB=ReadWriteIOUtils.readByte(inputStream)&0xFF;
        addCluster(L,((L<=R-interval_v)?(L+interval_v-1):R),cntB);
//        System.out.println("\t\t["+L+","+((L<=R-interval_v)?(L+interval_v-1):R)+"]"+":"+cntB);
//        debug_n+=cntB;
        L+=interval_v;
      }
      rankL=rankR;
      rankR = Math.min(rankL+stat_K,N-1);
      L=R;
      addSingleNum(R);
    }
//    System.out.println("\t\tcheck n: "+debug_n+" == "+N);
//    System.out.println("\t\tclusterLength "+clusterNum+"/"+clusterMinMaxSize.length/3);
  }

  public TDigestForPageStatReader(ByteBuffer byteBuffer) {
    this.N = ReadWriteIOUtils.readInt(byteBuffer);
    this.stat_K = ReadWriteIOUtils.readShort(byteBuffer);
    this.rnd = ReadWriteIOUtils.readShort(byteBuffer);// rnd=1..stat_K
//    System.out.println("\t\treading TDigest.\tN:"+N+"\tstat_K:"+stat_K+"\trnd:"+rnd);
    long L,R,interval_v;
    int rankL=0,rankR = Math.min(N,rnd)-1;
    int cntB;
    int debug_n = 0;
    L = ReadWriteIOUtils.readLong(byteBuffer);
    addSingleNum(L);
//    System.out.println("\t{"+L+"}");
//    debug_n++;
    if(rankR==0)
      rankR = Math.min(stat_K,N-1);
    while(rankL<rankR){
      R = ReadWriteIOUtils.readLong(byteBuffer);
//      System.out.println("\t{"+R+"}");
//      debug_n++;
      interval_v = L==R?1:(((R-L)>>>BUCKET_K)+1);
      while(L<=R){
        cntB=ReadWriteIOUtils.readByte(byteBuffer)&0xFF;
        addCluster(L,((L<=R-interval_v)?(L+interval_v-1):R),cntB);
//        System.out.println("\t\t["+L+","+((L<=R-interval_v)?(L+interval_v-1):R)+"]"+":"+cntB);
//        debug_n+=cntB;
        L+=interval_v;
      }
      rankL=rankR;
      rankR = Math.min(rankL+stat_K,N-1);
      L=R;
      addSingleNum(R);
    }
//    System.out.println("\t\tcheck n: "+debug_n+" == "+N);
//    System.out.println("\t\tclusterLength "+clusterNum+"/"+clusterMinMaxSize.length/3);
  }


}