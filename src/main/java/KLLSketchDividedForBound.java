import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;
import it.unimi.dsi.fastutil.longs.LongArrayList;

import java.util.*;

public class KLLSketchDividedForBound extends KLLSketchForQuantile {
  IntArrayList index;
  long[] allV,vInPage;
  LongArrayList v;
  int maxPageSizeOne;
  public long minV=Long.MAX_VALUE,maxV=Long.MIN_VALUE;

  public KLLSketchDividedForBound(int maxLevel, LongArrayList v) {
//    System.out.println("\t\t\tLongKLLSketchForDivide vSize:"+v.size());
    this.cntLevel = maxLevel + 1;
    N = v.size();
    levelPos = new int[cntLevel + 1];
    this.v = v;
  }

  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(N);
    sb.append((byte) cntLevel);
    sb.append((short) (levelPos[cntLevel] - levelPos[0]));
    for (int i = 0; i < cntLevel; i++) sb.append(levelPos[i]);
    return sb.toString();
  }

  private int getRankInPage(long v,int pageID,long[] VInPage,int L,int R){
    int rk=0;
    for(int i=maxPageSizeOne;i>0;i>>=1){
      if(L+rk+i<R&&VInPage[L+rk+i]<=v)rk+=i;
    }
    return rk+(VInPage[L+rk]<=v?1:0);
  }
  private int getDeltaRankInPage(long v,int lastRank,int wantWeight,long[] VInPage,int L,int R){
//    int leRk=0;
//    for(int i=maxPageSizeOne;i>0;i>>=1){
//      if(L+leRk+i<R&&VInPage[L+leRk+i]<v)leRk+=i;
//    }
//    leRk+=(VInPage[L+leRk]<=v?1:0);
    int leqRk=0;
    for(int i=maxPageSizeOne;i>0;i>>=1){
      if(L+leqRk+i<R&&VInPage[L+leqRk+i]<=v)leqRk+=i;
    }
    leqRk+=(VInPage[L+leqRk]<=v?1:0);
//    if(leqRk-lastRank>wantWeight)return Math.min(leRk-lastRank,wantWeight);
//    return Math.min(wantWeight,leqRk-lastRank);
    return leqRk-lastRank;
  }

  public List<LongKLLSketchForUnit> divideMemSketchByItemValue(IntArrayList pageStartIndex, LongArrayList pageMinV, LongArrayList pageMaxV){
//    System.out.println("\t\t\t???!!!");
//    allV = new long[v.size()];
//    v.toArray(allV);
    allV = v.toLongArray();
//    this.pageStartIndex = pageStartIndex;
    int pageNum = pageStartIndex.size();
    int[] pageIndexL = new int[pageNum],pageIndexR = new int[pageNum];
    for(int i=0;i<pageNum;i++) {
      int L = pageStartIndex.getInt(i), R = (i + 1 == pageNum ? (int) N : pageStartIndex.getInt(i + 1));
      Arrays.sort(allV, L, R);
      pageIndexL[i] = L;
      pageIndexR[i] = R;
      maxPageSizeOne = Math.max(maxPageSizeOne, Integer.highestOneBit(R-L));
    }
    vInPage = Arrays.copyOf(allV,(int)N);
    Arrays.sort(allV);
    minV = allV[0];
    maxV = allV[(int)N-1];
    index = new IntArrayList((int)N);
    for(int i=0;i<N;i++)index.add(i);
    fullyCompaction();

//    minV2Page = new Long2IntAVLTreeMap();
//    maxV2Page = new Long2IntAVLTreeMap(Comparator.comparingLong(x -> -x));

    IntHeapPriorityQueue minVinPage = new IntHeapPriorityQueue(pageNum,(x,y)->Long.compare(pageMinV.getLong(x),pageMinV.getLong(y)));
    PriorityQueue<Integer> maxVinPage = new PriorityQueue<>(pageNum, Comparator.comparingLong(pageMaxV::getLong));
    int[] lastItemRank = new int[pageNum];
//    ObjectHeapPriorityQueue<IntIntPair> dividingPage = new ObjectHeapPriorityQueue<>(pageNum,
//        (x,y)->(
//            x.getOne()==y.getOne()
//                ?Long.compare(pageMaxV.getLong(x.getTwo()),pageMaxV.getLong(y.getTwo()))
//                :Integer.compare(x.getOne(),y.getOne())
//            ));

    for(int i=0;i<pageNum;i++){
      minVinPage.enqueue(i);
    }
    int[] belongsTo = new int[levelPos[cntLevel]];
    int[] pageSketchSize = new int[pageNum], sameDeltaCount = new int[pageNum];
    int PAGE_SKETCH_SIZE = levelPos[cntLevel]/pageNum;
    int PAGE_ITEM_WEIGHT = 1<<(cntLevel-1);
//    System.out.println("\t\t??!! cntLevel:"+cntLevel+"\tPAGE_ITEM_WEIGHT:"+PAGE_ITEM_WEIGHT);
    List<LongKLLSketchForUnit> pageSketch = new ArrayList<>(pageNum);
    int itemLevel=0;
    for(int i=0;i<levelPos[cntLevel];i++){
      itemLevel+=(itemLevel<cntLevel-1&&i>=levelPos[itemLevel+1])?1:0;
      long itemValue = allV[index.getInt(i)];
      while(!minVinPage.isEmpty()&&pageMinV.getLong(minVinPage.firstInt())<=itemValue){
        maxVinPage.add(minVinPage.firstInt());
        minVinPage.dequeueInt();
      }
      while(!maxVinPage.isEmpty()&&pageMaxV.getLong(maxVinPage.peek())<itemValue){
        maxVinPage.poll();
      }
      int maxDeltaRank=Integer.MIN_VALUE,bestPage=-1;
      long bestDeltaV2MaxV=0;
      int sameDelta=0;
//      System.out.print(longToResult(itemValue)+"(");
      for(Iterator<Integer> it = maxVinPage.iterator();it.hasNext();){
        int cntPage = it.next();
//        int deltaRank = getRankInPage(itemValue,cntPage,vInPage,pageIndexL[cntPage],pageIndexR[cntPage]) - lastItemRank[cntPage];
        int deltaRank = getDeltaRankInPage(itemValue,lastItemRank[cntPage],PAGE_ITEM_WEIGHT,vInPage,pageIndexL[cntPage],pageIndexR[cntPage]);
        long cntDeltaV2MaxV = pageMaxV.getLong(cntPage)-itemValue;
        if(deltaRank>maxDeltaRank
            ||/*(deltaRank==maxDeltaRank&&pageSketchSize[cntPage]<pageSketchSize[bestPage])
            ||*/(deltaRank==maxDeltaRank&&/*pageSketchSize[cntPage]==pageSketchSize[bestPage]&&*/cntDeltaV2MaxV<bestDeltaV2MaxV)){
          sameDelta = deltaRank==maxDeltaRank?1:0;
          maxDeltaRank = deltaRank;
          bestPage = cntPage;
          bestDeltaV2MaxV = cntDeltaV2MaxV;
        }
//        System.out.print(deltaRank+",");
//        System.out.print("[Page"+cntPage+/*"index"+pageIndexL[cntPage]+"--"+pageIndexR[cntPage]+*/" Δ" + deltaRank+/*"  rk"+getRankInPage(itemValue,cntPage,vInPage,pageIndexL[cntPage],pageIndexR[cntPage])+*/"],");
      }
//      System.out.print(")");
//      lastItemRank[bestPage]+=maxDeltaRank;
      lastItemRank[bestPage]+=Math.min(PAGE_ITEM_WEIGHT, maxDeltaRank);
      belongsTo[i] = bestPage;
      pageSketchSize[bestPage]++;
      sameDeltaCount[bestPage]+=sameDelta;
//      System.out.print("_"+maxDeltaRank+"->"+bestPage+"\t\t");
//      if(pageSketchSize[bestPage]==PAGE_SKETCH_SIZE){
//        maxVinPage.remove(bestPage);
//      }
    }
//    System.out.print("\n\t\tItemLastRank:\t");
//    for(int i=0;i<pageNum;i++)System.out.print("\t\t"+lastItemRank[i]);System.out.println();
    for(int i=0;i<pageNum;i++)
      pageSketch.add(new LongKLLSketchForUnit(pageIndexR[i]-pageIndexL[i],cntLevel,pageSketchSize[i]));
    for(int i=0;i<cntLevel;i++){
      for(int j=0;j<pageNum;j++)pageSketch.get(j).levelPos[i+1]=pageSketch.get(j).levelPos[i];
      for(int j=levelPos[i];j<levelPos[i+1];j++)
        pageSketch.get(belongsTo[j]).addWhenDivide(i,allV[index.getInt(j)]);
    }
//    for(int i=0;i<pageNum;i++){
//      System.out.print("\t\t\t\trank In pageSketch:\t");
//      int lastRk=0;
//      double avgDelta=0;
//      for(int j=0;j<levelPos[cntLevel];j++)
//        if(belongsTo[j]==i) {
//          int cntRk = getRankInPage(allV[index.getInt(j)], i, vInPage, pageIndexL[i], pageIndexR[i]);
//          System.out.print("\t" + cntRk);
//          avgDelta+=1.0*(cntRk-lastRk)/pageSketchSize[i];
//          lastRk = cntRk;
//        }
//      System.out.println("\n\t\t\t\t\tAvgDelta="+avgDelta);
//    }
    return pageSketch;
  }

  public void sortIndex(){}


  public void fullyCompaction() {
    int cntPos = 0;
    int restN = (int)N;
//    System.out.println("\t\t"+ index.subList(cntPos,cntPos+restN));
    for (int i = 0; i < cntLevel-1; i++) {
      levelPos[i] = cntPos;
      if ((restN & 1) == 1) {
        /*int reserveP = cntPos + random.nextInt(restN);
        int reserveIndex = index.getInt(reserveP);
        for (int j = reserveP; j > cntPos; j--) index.set(j,index.getInt(j-1));
        index.set(cntPos,reserveIndex);*/ // reserve the first(min).
        cntPos++;
        restN--;
      }
      int odd = random.nextBoolean() ? 1 : 0;
      for (int j = 0; j < (restN >> 1); j++) {
//        odd = random.nextBoolean() ? 1 : 0;
        index.set(cntPos + j, index.getInt(cntPos + (j << 1 | (odd))));
      }
      restN >>= 1;
//      System.out.println("\t\t"+ index.subList(cntPos,cntPos+restN));
    }
    levelPos[cntLevel-1] = cntPos;
    levelPos[cntLevel] = cntPos + restN;
  }
}