import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactMultiQuantilesBestPrKLL {
    static int MULTI_QUANTILES;
    static boolean DEBUG_PRINT=false;
    int dataType=-233;
    static int startType=0,endType=0;
    static int pageN = 4096,batchN=4096;
    static int N = /*81000000*/8192*6713*2, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a;
    static DoubleArrayList prList=null,singlePrList=null;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    public XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom(233);

    static int DEBUG_COUNT_SIMULATION=0;

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextGaussian();
//                    i-0.5*Math.floor(i/pageN)*pageN;
//                    i;
            return;
        }
        if (dataType == 5) {
            for (int i = 0; i < N; i++)
                a[i] = random.nextGaussian();
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

//    public void prepareWorker(int maxSeriesByte) {
//        pageMinV = new double[pageNum];pageMaxV = new double[pageNum];
//
//        int enoughMemByte = pageN * 9;
//        for (int i = 0; i < pageNum; i++) {
//            pageMinV[i] = Long.MAX_VALUE;
//            pageMaxV[i] = Long.MIN_VALUE;
//            for (int j = 0; j < pageN; j++) {
//                double v =a[i*pageN+j];
//                pageMinV[i] = Math.min(pageMinV[i],v);
//                pageMaxV[i]=Math.max(pageMaxV[i],v);
//            }
////            worker.show();
//        }
//    }

    public int getValueActualRank(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(v<sortedA[mid])R=mid-1;
            else L=mid;
        }
        return L;
    }
    public int getValueLessThan(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(sortedA[mid]<v)L=mid;
            else R=mid-1;
        }
        return sortedA[L]<v?L:L-1;
    }
    public int getDeltaRank(double[]sortedA,int queryN, double v,int targetRank){
        int rank_L = getValueLessThan(sortedA,queryN,v)+1;
        int rank_R = getValueActualRank(sortedA,queryN,v);
//        System.out.println("\t\t\t"+targetRank+"\t\tresultLR:"+rank_L+"..."+rank_R+"\t\tresV:"+v);
        if(targetRank>=rank_L&&targetRank<=rank_R)return 0;
        else return targetRank<rank_L?(targetRank-rank_L):(targetRank-rank_R);
    }
    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    KLLSketchLazyEmptyForSimuCompact simuWorker;

    private double[] simulateIteration(double casePr,double lastPr,int depth,int n, int maxMemoryNum,int[] compactNum) {
        if (n <= 0) return new double[]{0,1.0};
        if (n <= maxMemoryNum) return
            new double[]{
//            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
                1.0,
                1.0};
        int maxERR = 0;
        for (int i = 0; i < compactNum.length; i++) maxERR += compactNum[i] << i;
        double bestSimuIter = 1e3, bestSimuPr = 1.0;

        double pr = bestSimuPr = lastPr;
        int prERR = KLLSketchLazyExact.queryRankErrBound(compactNum, pr);
        int successN=prERR * 2;
        int failN = (Math.min(n, maxERR * 2) - prERR) / 2;
//            KLLSketchLazyEmptyForSimuCompact simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*n, */maxMemoryNum);
        int[] succComNum,failComNum;
        if(failN<successN){
            failComNum = simuWorker.simulateCompactNumGivenN(failN);
            succComNum = simuWorker.simulateCompactNumGivenN(successN);
        }else{
            succComNum = simuWorker.simulateCompactNumGivenN(successN);
            failComNum = simuWorker.simulateCompactNumGivenN(failN);
        }
        bestSimuIter = 1 + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum,succComNum)[0] + (1 - pr) * (1 + simulateIteration(casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum,failComNum)[0]);

        DEBUG_COUNT_SIMULATION++;
        return new double[]{bestSimuIter,bestSimuPr};
    }
//    private double[] evaluatePr(int maxMemoryNum,int prID,int succN,int failN,int[] succComNum,int[]failComNum){
//        double totPr=prList.getDouble(prID),singlePr=singlePrList.getDouble(prID);
//        double[] simulateResult= new double[3];
//        double[] successResult = simulateIteration(totPr,Pr,1,succN,maxMemoryNum,succComNum);
//        double[] failResult = simulateIteration(0*(1-Pr),Pr,1,failN,maxMemoryNum,failComNum);
//        simulateResult[0] = Pr*successResult[0]+(1-Pr)*(1+failResult[0]);
//        simulateResult[1] = successResult[1];
//        simulateResult[2] = failResult[1];
////        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\t\t\tcntPR:"+Pr+"\tsuccessN:\t"+succN+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateResult[0]+"\tsimu_nextSuccessPr:"+simulateResult[1]);
//        return simulateResult;
//    }
    private int findBestPrID(KLLSketchLazyExact sketch, int maxMemoryNum,long rk1,long rk2,double[] deterministic_result,int MULTI_QUANTILES){
        double bestEstiNum=1e9,nextSuccessPr=1.0;
        int bestPrId = 0;
        int[] successN = new int[prList.size()],failN = new int[prList.size()];
        int[][] succComNum= new int[prList.size()][],failComNum=new int[prList.size()][];
        for(int i=0;i<prList.size();i++){
            double totPr=prList.getDouble(i),singlePr=singlePrList.getDouble(i);
            double[] cntResult = sketch.findResultRange(rk1,rk2,singlePr);
            int rkValL = (int)cntResult[2],rkValR=(int)cntResult[3],prErrL=(int)cntResult[4],prErrR=(int)cntResult[5];
            int tmpSuccessN = Math.max(rkValR-rkValL,prErrL+prErrR);
            tmpSuccessN+=(prErrL+prErrR)/16;
            if(tmpSuccessN<=maxMemoryNum)
                tmpSuccessN+=(prErrL+prErrR)/16;
            int tmpFailN = (((int)deterministic_result[3]-(int)deterministic_result[2])-tmpSuccessN)/2;
            successN[i]=tmpSuccessN;
            failN[i]=Math.max(tmpSuccessN,tmpFailN);
        }
        //KLLSketchLazyEmptyForSimuCompact
        simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*(int)sketch.getN()/2, */maxMemoryNum);
        for(int i=0;i<prList.size();i++)
            succComNum[i]=simuWorker.simulateCompactNumGivenN(successN[i]);
        for(int i=prList.size()-1;i>=0;i--)
            failComNum[i]=simuWorker.simulateCompactNumGivenN(failN[i]);

        for(int i=0;i<prList.size();i++) {
            double cntPrResult = evaluatePrForMultiQuantiles(maxMemoryNum, i, successN[i],failN[i],succComNum[i],failComNum[i],MULTI_QUANTILES);
            if(cntPrResult<=bestEstiNum){
                bestEstiNum = cntPrResult;
                bestPrId = i;
            }
//            System.out.println("\t\t\t\t\tcntPR:"+prList.getDouble(i)+"\tsuccessN:\t"+successN[i]+"\t\tfailN:\t"+failN[i]+"\t\testiIter:\t"+bestEstiNum);
        }
        if(DEBUG_PRINT)System.out.println("\tbestTotPr:\t"+prList.getDouble(bestPrId)+"\t\tbestSinglePr:\t"+singlePrList.getDouble(bestPrId)+"\t\testiIter:"+bestEstiNum+"\tnextSuccessN:\t"+successN[bestPrId]);
        return bestPrId;
    }


    private double simulateIterationForMultiQuantiles(double casePr,int prID,int depth,int n, int maxMemoryNum,int[] compactNum) {
        if (n <= 0) return 0;
        if (n <= maxMemoryNum) return
//            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
                1.0;
        int maxERR = 0;
        for (int i = 0; i < compactNum.length; i++) maxERR += compactNum[i] << i;
        double bestSimuIter = 1e3;

        double totPr=prList.getDouble(prID),singlePr=singlePrList.getDouble(prID);
        int prERR = KLLSketchLazyExact.queryRankErrBound(compactNum, singlePr);
        int successN=prERR * 2;
        int failN = (Math.min(n, maxERR * 2) - prERR) / 2;
        failN=Math.max(failN,successN);

//            KLLSketchLazyEmptyForSimuCompact simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*n, */maxMemoryNum);
        int[] succComNum,failComNum;
        succComNum = simuWorker.simulateCompactNumGivenN(successN);
        failComNum = simuWorker.simulateCompactNumGivenN(failN);
        bestSimuIter = 1 + totPr * simulateIterationForMultiQuantiles(casePr * totPr, prID, depth + 1, successN, maxMemoryNum,succComNum) + (1 - totPr) * (1 + simulateIterationForMultiQuantiles(casePr * (1 - totPr), prID, depth + 1, failN, maxMemoryNum,failComNum));

        DEBUG_COUNT_SIMULATION++;
        return bestSimuIter;
    }
    private double evaluatePrForMultiQuantiles(int maxMemoryNum,int prID,int succN,int failN,int[] succComNum,int[]failComNum,int MULTI_QUANTILES){
        double totPr=prList.getDouble(prID),singlePr=singlePrList.getDouble(prID);
        double successResult = simulateIterationForMultiQuantiles(totPr,prID,1,succN,maxMemoryNum,succComNum);
        double failResult = simulateIterationForMultiQuantiles(0*(1-totPr),prID,1,failN,maxMemoryNum,failComNum);
        double simulateResult = totPr*successResult+(1-totPr)*(1+failResult);
//        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\t\t\tcntPR:"+Pr+"\tsuccessN:\t"+succN+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateResult[0]+"\tsimu_nextSuccessPr:"+simulateResult[1]);
        return simulateResult;
    }
    private int findBestPrIDForMultiQuantiles(KLLSketchLazyExact sketch, int maxMemoryNum,int MULTI_QUANTILES){
        for(int i=0;i<prList.size();i++)
            singlePrList.set(i,Math.pow(prList.getDouble(i),1.0/MULTI_QUANTILES));
        double bestEstiNum=1e9;
        int bestPrId = 0,prNum=prList.size();
        int[] successN = new int[prNum],failN = new int[prNum];
        int[][] succComNum= new int[prNum][],failComNum=new int[prNum][];
        int n=(int)sketch.getN(),maxERR=sketch.getMaxErr();
        for(int i=0;i<prNum;i++){
            int prERR = sketch.queryRankErrBound(singlePrList.getDouble(i));
            successN[i]=prERR * 2;
            failN[i]=(Math.min(n, maxERR * 2) - prERR) / 2;
            failN[i]=Math.max(failN[i],successN[i]);
        }
        //KLLSketchLazyEmptyForSimuCompact
        simuWorker = new KLLSketchLazyEmptyForSimuCompact(maxMemoryNum);
        for(int i=0;i<prNum;i++)
            succComNum[i]=simuWorker.simulateCompactNumGivenN(successN[i]);
        for(int i=prNum-1;i>=0;i--)
            failComNum[i]=simuWorker.simulateCompactNumGivenN(failN[i]);

        for(int i=0;i<prNum;i++) {
            double cntPrResult = evaluatePrForMultiQuantiles(maxMemoryNum, i, successN[i],failN[i],succComNum[i],failComNum[i],MULTI_QUANTILES);
            if(cntPrResult<=bestEstiNum){
                bestEstiNum = cntPrResult;
                bestPrId = i;
            }
//            System.out.println("\t\t\t\t\tcntPR:"+prList.getDouble(i)+"\tsuccessN:\t"+successN[i]+"\t\tfailN:\t"+failN[i]+"\t\tsimu_iter:\t"+cntNum);
        }
        if(DEBUG_PRINT)System.out.println("bestTotPr:\t"+prList.getDouble(bestPrId)+"\t\tbestSinglePr:\t"+singlePrList.getDouble(bestPrId)+"\t\testiIter:"+bestEstiNum+"\tnextSuccessN:\t"+successN[bestPrId]);
        return bestPrId;
    }





    public void testError(int queryN, int maxMemoryByte,int MULTI_QUANTILES) throws IOException {
        if (prList == null) {
            prList = new DoubleArrayList();
            for (double tmp = 0.70; tmp < 0.9 - 1e-6; tmp += 0.02) prList.add(tmp);
            for (double tmp = 0.90; tmp < 0.99 - 1e-6; tmp += 0.01) prList.add(tmp);
            prList.add(0.99);
            prList.add(0.995);
            singlePrList = new DoubleArrayList(prList);
        }
        this.batchN=MULTI_QUANTILES;
//        this.MULTI_QUANTILES=MULTI_QUANTILES;
        double[] query_q = new double[MULTI_QUANTILES];
        final double query_delta_q = 1.0 / (MULTI_QUANTILES+1);
        for (int i = 0; i < MULTI_QUANTILES; i++) query_q[i] = query_delta_q * (i + 1);
        long full_time = 0, merge_page_time = 0;
        double avg_iteration = 0;
//        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum = 0, ALLPrCount = 0, FailCount = 0, IterCount = 0;

//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        for (int T = 0; T < TEST_CASE; T++) {
//            if(T%100==0) System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L = LL[T], R = RR[T];

//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);

            double ratioPerQuery = 1.0 / TEST_CASE;
            int[] query_rank1=new int[MULTI_QUANTILES],query_rank2=new int[MULTI_QUANTILES];
            double[][] deterministic_result=new double[MULTI_QUANTILES][], iterate_result=new double[MULTI_QUANTILES][];
            for(int qid=0;qid<MULTI_QUANTILES;qid++) {
                query_rank1[qid] = (int) Math.floor(query_q[qid] * (queryN - 1) + 1);
                query_rank2[qid] = (int) Math.ceil(query_q[qid] * (queryN - 1) + 1);
                deterministic_result[qid]=new double[]{-Double.MAX_VALUE, Double.MAX_VALUE};
                iterate_result[qid]=new double[]{-Double.MAX_VALUE, Double.MAX_VALUE};
            }
            int MMP = 0;
            while (true) {
                boolean calc_finished=true;
                for(int qid=0;qid<MULTI_QUANTILES;qid++)
                    if(deterministic_result[qid][0] < deterministic_result[qid][1] && deterministic_result[qid].length != 3)
                        calc_finished=false;
                if(calc_finished)break;
                if(++MMP>15){
                    System.out.println("\t\t\t\t\t\titerate fail.\titer:"+MMP);
                    for(int qid=0;qid<MULTI_QUANTILES;qid++)
                        if(deterministic_result[qid][0] < deterministic_result[qid][1] && deterministic_result[qid].length != 3)
                            System.out.println("\t\t\t\t\t\t\t\t\t\tremaining_q:\t"+query_q[qid]);
                    return;
                }
                    System.out.println("\n\tnew iter.??\t"+iterate_result[0][0]+"..."+iterate_result[0][1]);
                avg_iteration += ratioPerQuery;
                IterCount += 1;

                if(MULTI_QUANTILES>1&&MMP==1){// first iter.
                    KLLSketchLazyExactPriori cntWorker = new KLLSketchLazyExactPriori(maxMemoryByte);
//                    FastKLLSketchLazyForBound cntWorker = new FastKLLSketchLazyForBound(maxMemoryByte);
                    int CountOfLessThanValL = 0;
                    for (int i = L; i < R; i++) cntWorker.update(dataToLong(a[i]));
                    if (DEBUG_PRINT)
                        System.out.println("\t\t\t\t\t\titerate success." + "\t\tcntN:" + cntWorker.getN() + "\t\tcntIter:" + MMP);
                    for(int qid=0;qid<MULTI_QUANTILES;qid++)
                        deterministic_result[qid] = cntWorker.getFilter(0,0,0,0,query_rank1[qid], query_rank2[qid], 1.0);

                    int[] compactNum =cntWorker.getRelatedCompactNum();
                    PrioriBestPrHelper helper = new PrioriBestPrHelper(maxMemoryByte/8,queryN,compactNum,0,MULTI_QUANTILES);
                    double bestPr = helper.findBestPr(5e-4,5e-4,5e-1,queryN)[0];
                    ALLPrSum += bestPr;
                    ALLPrCount += 1;
                    for(int qid=0;qid<MULTI_QUANTILES;qid++)
                        iterate_result[qid] = cntWorker.getFilter(0,0,0,0,query_rank1[qid], query_rank2[qid], Math.pow(bestPr,1.0/MULTI_QUANTILES));
                    System.out.println("\t\t\tbestPr:\t"+bestPr);
                }else {
                    int rest_multi_quantiles = 0;
                    int[] rest_query_id=new int[MULTI_QUANTILES];
                    int[] rest_query_rank1=new int[MULTI_QUANTILES],rest_query_rank2=new int[MULTI_QUANTILES];
                    double[][] rest_deterministic_result=new double[MULTI_QUANTILES][], rest_iterate_result=new double[MULTI_QUANTILES][];
                    for(int qid=0;qid<MULTI_QUANTILES;qid++)
                        if(deterministic_result[qid][0] < deterministic_result[qid][1] && deterministic_result[qid].length != 3) {
                            rest_query_rank1[rest_multi_quantiles] = query_rank1[qid];
                            rest_query_rank2[rest_multi_quantiles] = query_rank2[qid];
                            rest_deterministic_result[rest_multi_quantiles] = deterministic_result[qid];
                            rest_iterate_result[rest_multi_quantiles] = iterate_result[qid];
                            rest_query_id[rest_multi_quantiles]=qid;
                            rest_multi_quantiles++;
                        }

                    KLLSketchLazyExactPriori[] cntWorker = new KLLSketchLazyExactPriori[rest_multi_quantiles];

                    double minValL = rest_iterate_result[0][0], maxValR = rest_iterate_result[0][1];
                    for(int i=0;i<rest_multi_quantiles;i++){
                        cntWorker[i]=new KLLSketchLazyExactPriori(maxMemoryByte/rest_multi_quantiles);
                        minValL=Math.min(minValL,rest_iterate_result[i][0]);
                        maxValR=Math.max(maxValR,rest_iterate_result[i][1]);
                    }
//                    if(DEBUG_PRINT)
                        System.out.println("\n\n\t\t\trest_multi_quantiles:\t"+rest_multi_quantiles);

                    int[] CountOfLessThanValL = new int[rest_multi_quantiles+1];
                    int[] CountOfValL = new int[rest_multi_quantiles+1];
                    int[] CountOfValR = new int[rest_multi_quantiles+1];
                    long[] valBuf=new long[batchN];
                    int bufSize=0;
                    for (int dataI = L; dataI < R; dataI++) {
                        valBuf[bufSize++]=dataToLong(a[dataI]);
                        if(bufSize==batchN) {
                            Arrays.sort(valBuf);
                            int index=0;
                            for(long dataL:valBuf){
                                double dataD=longToResult(dataL);
                                while(index<rest_multi_quantiles&&dataD>(rest_iterate_result[index][1]))index++;
                                if(index==rest_multi_quantiles)break;
                                int i=index;
                                for(;i<rest_multi_quantiles&&dataD>=(rest_iterate_result[i][0]);i++)
                                    if(dataD<(rest_iterate_result[i][0])) {
                                        CountOfLessThanValL[i]++;
                                        CountOfLessThanValL[i+1]--;
                                    }else
                                    if((rest_iterate_result[i][0])<=dataD&&dataD<=(rest_iterate_result[i][1])) {
                                        if(dataD==rest_iterate_result[i][0])CountOfValL[i]++;
                                        else if(dataD==rest_iterate_result[i][1])CountOfValR[i]++;
                                        else cntWorker[i].update(dataL);
                                    }
                                CountOfLessThanValL[i]++;
                            }
                            bufSize=0;
                        }
                    }
                    Arrays.sort(valBuf,0,bufSize);
                    int index=0;
                    for(int bufID=0;bufID<bufSize;bufID++){
                        long dataL=valBuf[bufID];
                        double dataD=longToResult(dataL);
                        while(index<rest_multi_quantiles&&dataD>(rest_iterate_result[index][1]))index++;
                        if(index==rest_multi_quantiles)break;
                        int i=index;
                        for(;i<rest_multi_quantiles&&dataD>=(rest_iterate_result[i][0]);i++)
                            if(dataD<(rest_iterate_result[i][0])) {
                                CountOfLessThanValL[i]++;
                                CountOfLessThanValL[i+1]--;
                            }else
                            if((rest_iterate_result[i][0])<=dataD&&dataD<=(rest_iterate_result[i][1])) {
                                if(dataD==rest_iterate_result[i][0])CountOfValL[i]++;
                                else if(dataD==rest_iterate_result[i][1])CountOfValR[i]++;
                                else cntWorker[i].update(dataL);
                            }
                        CountOfLessThanValL[i]++;
                    }
                    bufSize=0;
                    for(int i=1;i<rest_multi_quantiles;i++)
                        CountOfLessThanValL[i] += CountOfLessThanValL[i - 1];

                    IntArrayList toEstimateID = new IntArrayList(rest_multi_quantiles);
                    for(int i=0;i<rest_multi_quantiles;i++) {
//                        if(DEBUG_PRINT)cntWorker[i].show();
                        int qid=rest_query_id[i];
                        int cntRank1 = rest_query_rank1[i] - CountOfLessThanValL[i];
                        int cntRank2 = rest_query_rank2[i] - CountOfLessThanValL[i];
//                    System.out.println("\t\t\t\t\t\tcntRank:"+cntRank1+" "+cntRank2);
                        if (cntRank1 <= 0 || cntRank2 > CountOfValL[i]+CountOfValR[i]+cntWorker[i].getN()) { // iteration failed.
                            if (DEBUG_PRINT) System.out.println("\t\t\t\t\t\titerate fail." + "\t\tcntIter:" + MMP+"\t\tcntQ:"+query_q[qid]);
                            if (cntRank1 <= 0)
                                iterate_result[qid] = new double[]{rest_deterministic_result[i][0], rest_iterate_result[i][0]};
                            else iterate_result[qid] = new double[]{rest_iterate_result[i][1], rest_deterministic_result[i][1]};
                            deterministic_result[qid] = iterate_result[qid];
                            if (deterministic_result[qid][0] == deterministic_result[qid][1]) continue;
                            FailCount += 1;
//                        IterCount-=1;
                            continue;
                        }
                        if (DEBUG_PRINT)
                            System.out.println("\t\t\t\t\t\titerate success." + "\t\tcntN:" + cntWorker[i].getN() + "\t\tcntIter:" + MMP);
                        deterministic_result[qid] = cntWorker[i].getFilter(CountOfValL[i],CountOfValR[i],iterate_result[i][0],iterate_result[i][1],cntRank1, cntRank2, 1.0);
                        if (deterministic_result[qid].length == 3 || deterministic_result[qid][0]==deterministic_result[qid][1]) {
                            iterate_result[qid] = deterministic_result[qid];
                            continue;
                        }
                        toEstimateID.add(i);
                    }
                    if (toEstimateID.size() > 0) {
                        int worstID = toEstimateID.getInt(0);
                        for(int i:toEstimateID)
                            if(cntWorker[i].getN()>cntWorker[worstID].getN())worstID=i;
                        int[] compactNum =cntWorker[worstID].getRelatedCompactNum();
                        PrioriBestPrHelper helper = new PrioriBestPrHelper(maxMemoryByte/8,(int)cntWorker[worstID].getN(),compactNum,0,toEstimateID.size());
                        double bestPr = helper.findBestPr(5e-4,5e-4,5e-1,(int)cntWorker[worstID].getN())[0];
                        System.out.println("\t\t\tbestPr:\t"+bestPr);
                        System.out.println("\t\t\tworst com:\t"+ Arrays.toString(compactNum) +"\t\tworstN:\t"+(int)cntWorker[worstID].getN());

                        for(int i:toEstimateID) {
                            int qid=rest_query_id[i];
                            iterate_result[qid] = cntWorker[i].getFilter(CountOfValL[i],CountOfValR[i],iterate_result[i][0],iterate_result[i][1],
                                rest_query_rank1[i]-CountOfLessThanValL[i], rest_query_rank2[i]-CountOfLessThanValL[i], Math.pow(bestPr,1.0/toEstimateID.size()));
//                            deterministic_result[qid] = cntWorker[i].findResultRange(rest_query_rank1[i]-CountOfLessThanValL[i], rest_query_rank2[i]-CountOfLessThanValL[i], 1.0);
                        }

                        ALLPrSum += bestPr;
                        ALLPrCount += 1;
                    }
//                    System.out.println("\t\t\t\tbestPr:"+bestPr+"\t\t\tcntN:"+cntWorker.getN());
//                    System.out.println("\t\t\t\t\tcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);
                }
            }

//            double[] query_a = new double[queryN];
//            System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);
//            for(int qid=0;qid<MULTI_QUANTILES;qid++) {
//                double exact_calc_v = (deterministic_result[qid][0] + deterministic_result[qid][1]) * 0.5;
//                double exact_answer_v = (query_a[query_rank1[qid]-1]+query_a[query_rank2[qid]-1]) * 0.5;
////                if (DEBUG_PRINT)
////                    System.out.println("\tFINISH CALC\t\texact_calc_v=:\t" + exact_calc_v +"\tanswer_v=:\t"+exact_answer_v+ "\tqueryQ:\t"+query_q[qid]);
//                    if(Math.abs(exact_answer_v-exact_calc_v)/Math.max(Math.abs(exact_answer_v),Math.abs(exact_calc_v))>1e-9)
//                        System.out.println("\t\t!! CALC ERROR\n");
////                if(Math.abs(exact_quantile_v-N*query_q[qid])>1.0){
////                    System.out.println("!!!!!!!!!!ERROR!!!!");
////                    return;
////                }
//            }
        }

//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("MultiQuantiles!\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + maxMemoryByte + "\tMultiQuantiles:\t"+MULTI_QUANTILES+"\t\tavg_iteration:\t" + avg_iteration + "\t\tavgPrChose:" + ALLPrSum / ALLPrCount + "\t\tavgFailRate:" + (FailCount / IterCount) + "\tusedNorDistri:\t" + KLLSketchLazyExact.DEBUG_COUNT_USENORMAL + "\tnot_cached_calc_pr_err:" + KLLSketchLazyExact.DEBUG_COUNT_USEPRERR + "\t\tSIMU_COUNT:" + DEBUG_COUNT_SIMULATION);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactMultiQuantilesBestPrKLL main;

        main = new MainForExactMultiQuantilesBestPrKLL();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
//            for (int queryN : new int[]{/*1000000,1500000,2000000,2500000,3000000,*/3500000,4000000,4500000,5000000/*,5500000,6000000/**/})
            for (int MULTI_QUANTILES : new int[]{/*1,1,*/500/*,4,6,10,20,40,60,100,200,400,600,1000,2000,4000,6000,10000/**/})
            for (int queryN : new int[]{100000000/*,20000000,30000000,40000000,50000000,60000000,70000000,80000000*/})
                for (int query_mem : new int[]{1024*1024*8})
//                    for (int query_mem : new int[]{1024*128})
                    main.testError(queryN, query_mem,MULTI_QUANTILES);
                main.RESULT_LINE++;
//            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
