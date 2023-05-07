import org.apache.iotdb.tsfile.utils.ReadWriteIOUtils;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Random;

// based on KLL Sketch in DataSketch. See https://github.com/apache/datasketches-java/tree/master/src/main/java/org/apache/datasketches/kll
// This is an implementation for long data type.
// We changed the behaviour of serialization to reduce the serialization size. In our situation, there will be few update operations after deserialization.
public class OldLongKLLSketch {
  int maxN,N,maxSerializeNum;
  int K;
  long[] num;
  boolean level0Sorted;
  int maxLevel;
  int[] levelPos,levelMaxSize;
  long XORSHIFT=new Random().nextInt();//0x2333333319260817L;
//  Random test_random = new Random();


  public OldLongKLLSketch(int maxN, int levelNUM, int maxMemoryByte, int maxSerializeByte) {
    int maxTotalSize=0;
    this.maxN = maxN;
    this.N=0;
    maxLevel=levelNUM;//calcLevelNum(maxN,maxSerializationByte);
    this.maxSerializeNum = (maxSerializeByte-9-maxLevel)/8;
    K =(maxN-1)/(1<<(maxLevel-1))+1;
    levelPos = new int[maxLevel+1];
    levelMaxSize = new int[maxLevel];
    for(int i=maxTotalSize=0;i<maxLevel;i++){
      levelMaxSize[i] = Math.max(8,(int)Math.round((K*Math.pow(2.0/3,maxLevel-i-1))));
      maxTotalSize += levelMaxSize[i];
//      System.out.println("\t\t\tLevelMaxSize:"+levelMaxSize[i]);
    }
    int arraySize = maxTotalSize/*+7000*/;//maxTotalSize;//Math.max(maxTotalSize,K*3);
    arraySize=(maxMemoryByte/8-maxLevel*2-5);
    arraySize = Math.min(arraySize, 1<<22);
    int extraSize = arraySize-maxTotalSize;
//    {
//      for (int i = 0; i < maxLevel - 1; i++)
//        levelMaxSize[i] = (int) (levelMaxSize[i] * (1.0 + 1.0 * extraSize / maxTotalSize));
//    }
//    {
//      for (int i = maxLevel-1-1; i >=0; i--) {
//        int delta=Math.min(extraSize/(i+1),(K<<(maxLevel-i-1))-levelMaxSize[i]);
//        levelMaxSize[i]+=delta;
//        extraSize-=delta;
//      }
//    }
//    for(int sameLevel=1,remaining=extraSize;sameLevel<maxLevel;sameLevel++){
//      {
//        int need = 0;
//        for (int i = 0; i < sameLevel; i++) need += levelMaxSize[sameLevel] - levelMaxSize[i];
//        if (need <= remaining) {
//          remaining -= need;
//          for (int i = 0; i < sameLevel; i++) levelMaxSize[i] = levelMaxSize[sameLevel];
//        } else {
//          for (int i = 0; i < sameLevel; i++) levelMaxSize[i] += remaining / sameLevel;
//          break;
//        }
//      }
//    }
    {
      maxTotalSize = 0;
      int newK = K;
      for (int addK = 1 << 25; addK > 0; addK >>>= 1) {
        int need = 0;
        for (int i = 0; i < maxLevel - 1; i++)
          need += Math.min(K<<(maxLevel-i-1), Math.max(8, (int) Math.round(((newK + addK) * Math.pow(2.0 / 3, maxLevel - i - 1))))) - levelMaxSize[i];
        if (need <= extraSize) newK += addK;
      }
      for (int i = 0; i < maxLevel - 1; i++) {
        levelMaxSize[i] = Math.min(K << (maxLevel - i - 1), Math.max(8, (int) Math.round((newK * Math.pow(2.0 / 3, maxLevel - i - 1)))));
        maxTotalSize += levelMaxSize[i];
      }
//      System.out.println("\t\t\t\t\t"+maxTotalSize+"???"+maxMemoryByte);
      arraySize=Math.min(maxN,maxTotalSize);
    }

    num = new long[arraySize];
    for(int i=0;i<=maxLevel;i++)levelPos[i] = arraySize;
    level0Sorted = false;
//    System.out.println("\t\tArraySize:"+arraySize);
  }


  private int getMaxUnEmptyLevel(){
    for(int i=maxLevel-1;i>=0;i--)if(levelPos[i]!=levelPos[i+1])return i;
    return 0;
  }

  public String toString(){
    final StringBuilder sb = new StringBuilder();
    sb.append(N);
    sb.append(maxN);
    sb.append((byte)maxLevel);
    sb.append((short)(levelPos[maxLevel]-levelPos[0]));
    for(int i=0;i<maxLevel;i++)sb.append(levelPos[i]);
    return "";
  }

  public int serialize(OutputStream outputStream) throws IOException {
    compactBeforeSerialization();
    int byteLen = 0;
    byteLen+=ReadWriteIOUtils.write(N, outputStream);
    byteLen+=ReadWriteIOUtils.write(maxN, outputStream);
    byteLen+=ReadWriteIOUtils.write((byte)maxLevel, outputStream);
    if(K<100)
      for(int i=0;i<maxLevel;i++)byteLen+=ReadWriteIOUtils.write((byte)(levelPos[i+1]-levelPos[i]), outputStream);
    else
      for(int i=0;i<maxLevel;i++)byteLen+=ReadWriteIOUtils.write((short)(levelPos[i+1]-levelPos[i]), outputStream);
    for(int i=levelPos[0];i<levelPos[maxLevel];i++)ReadWriteIOUtils.write(num[i],outputStream);
    return byteLen;
  }
  private int calcLevelNum(int maxN, int maxByte){
//    return 3;
    for(int level=5;level<=20;level++){
      int numInTopLevel = (maxN>>>(level-1));
      int otherLevelNotEmpty = level-1;// Long.bitCount(actualN&((1<<(level-1))-1));
      int tmpByte = 4+8*2+1+otherLevelNotEmpty*8+numInTopLevel*8; // actualN,minmax,level, otherLevel, TopLevel
      System.out.println("\t\t\t\t "+level+"   byte:"+tmpByte);
      if(tmpByte<=maxByte)
        return level;
    }
    return -1;
  }
  public void update(long x){ // signed long
    if(levelPos[0]==0)compact();
    num[--levelPos[0]] = x;
    N++;
    level0Sorted = false;
//    if(N>=maxN*0.9)
//      compact();
//      compactBeforeSerialization();
//    if(N<=10)System.out.println("\t\t\t"+x);
  }
  public void show(){
    for(int i=0;i<maxLevel;i++){
      System.out.print("\t");
      System.out.print("["+(levelPos[i+1]-levelPos[i])+"]");
      System.out.print("\t");
    }System.out.println();
    for(int i=0;i<maxLevel;i++){
      System.out.print("\t");
      for(int j=levelPos[i];j<levelPos[i+1];j++)
        System.out.print(num[j]+",");
      System.out.print("\t");
    }System.out.println();
    System.out.println("\t\tCOMPACT_SIZE:"+(4+8*2+1+(maxLevel-1)*8+(int)((maxN>>>(maxLevel-1)))*8)+"\t//maxMemNum:"+levelPos[maxLevel]+"maxSeriNum:"+maxSerializeNum+"\t//N:"+N);
//    System.out.println("\t\t\t01:"+ZEROONESUM+"/"+ZEROONENUM+":  "+(1.0*ZEROONENUM/ZEROONESUM));

    for(int i=0;i<maxLevel;i++)System.out.print("\t\t"+levelMaxSize[i]+"\t");System.out.println();
    System.out.println("-------------------------------------------------------");
  }




  private int getNextRand01(){ // xor shift *
    XORSHIFT^=XORSHIFT>>>12;
    XORSHIFT^=XORSHIFT<<25;
    XORSHIFT^=XORSHIFT>>>27;
    return (int) ((XORSHIFT*0x2545F4914F6CDD1DL)&1);
//    return test_random.nextInt()&1;
  }
//  long ZEROONENUM=0,ZEROONESUM=0;
  private void randomlyHalveDownToLeft(int L,int R){
//    System.out.println("\t\t\t\t01:"+getNextRand01());
//    ZEROONENUM+=getNextRand01();
//    ZEROONESUM++;
    int delta = getNextRand01();
    int mid = (L+R)>>>1;
    for(int i=L,j=L;i<mid;i++,j+=2)
      num[i]=num[j+delta];
  }
  private void mergeSortWithSpace(int L1,int mid,int L2,int R2){
    int p1 = L1, p2=L2, cntPos=mid;
    while(p1<mid||p2<R2){
      if(p1<mid&&(p2==R2||num[p1]<num[p2]))num[cntPos++]=num[p1++];
      else num[cntPos++]=num[p2++];
    }
  }

  private void compactOneLevel(int level,int numToReduce){ // compact half of data when numToReduce is small
    if(numToReduce<=0)return;
    int L1 = levelPos[level], R1 = levelPos[level+1]; // [L,R)
    if(level==0&&!level0Sorted) {
      Arrays.sort(num, L1, R1);
      level0Sorted = true;
    }
    L1 += (R1-L1)&1;
    if(L1==R1)return;
//    boolean FLAG = false;
    if((R1-L1)/2>numToReduce){
      L1=R1-numToReduce*2;
//      FLAG=true;
//      System.out.println("----------------------------");show();
    }

    randomlyHalveDownToLeft(L1,R1);

    int mid=(L1+R1)>>>1;
    mergeSortWithSpace(L1,mid,levelPos[level+1],levelPos[level+2]);
    levelPos[level+1]=mid;
    int newP = levelPos[level+1]-1,oldP=L1-1;
    for(int i=oldP;i>=levelPos[0];i--)
      num[newP--]=num[oldP--];

    levelPos[level]=levelPos[level+1]-(L1-levelPos[level]);
    int numReduced=(R1-L1)>>>1;
    for(int i=level-1;i>=0;i--)levelPos[i]+=numReduced;
//    if(FLAG){
//      show();
//      System.out.println("----------------------------");
//    }
  }

  public void compact(){
    int compactLevel=maxLevel-1;
    for(int i=0;i<maxLevel;i++)
    if(levelPos[i+1]-levelPos[i] > levelMaxSize[i]){
      compactLevel = i;
      break;
    }
    if(compactLevel>=maxLevel-1){
//      System.out.println("\t\t[MYKLL ERROR] compactLevel wrong!");
      return;
    }
    compactOneLevel(compactLevel,Integer.MAX_VALUE);
  }
  public void compactBeforeSerialization(){
//    int topLevel = getMaxUnEmptyLevel();
    for (int i = 0; i < maxLevel - 1; i++) {
//      System.out.println("\t\t\t\t"+i);
      if (levelPos[i + 1] - levelPos[i] >= 2)
        compactOneLevel(i,levelPos[maxLevel]-levelPos[0] - maxSerializeNum);
      if(levelPos[maxLevel]-levelPos[0]<=maxSerializeNum)
        break;
    }
  }

  public int findRankInLevel(int level,long v){
    int L = levelPos[level],R = levelPos[level+1];
    if(level==0&&!level0Sorted) {
      Arrays.sort(num, L, R);
      level0Sorted = true;
    }
    R--;
    if(L>R || num[L]>=v)return 0;
    while(L<R){
      int mid=(L+R+1)>>1;
      if(num[mid]<v)L=mid;else R=mid-1;
    }return (L-levelPos[level]+1)*(1<<level);
  }


  public int getApproxRank(long v){
    int approxRank = 0;
    for(int i=0;i<maxLevel;i++) {
      approxRank+=findRankInLevel(i,v);
//      for (int j = levelPos[i]; j < levelPos[i + 1]; j++)
//        if (num[j] < v) approxRank += 1 << i;
    }
    return approxRank;
  }

//  public int getLowerBound(int rank){
//    compactBeforeSerialization();
//
//  }


}