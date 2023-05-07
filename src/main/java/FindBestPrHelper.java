import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

public class FindBestPrHelper {
  static KLLSketchLazyEmptyForSimuCompact simuWorker;
  boolean DEBUG_PRINT=false;
  static DoubleArrayList prList;
  static{
    prList = new DoubleArrayList();
    for(double tmp=0.70;tmp<0.9-1e-6;tmp+=0.02)prList.add(tmp);
    for(double tmp=0.90;tmp<0.99-1e-6;tmp+=0.01)prList.add(tmp);
    prList.add(0.99);
    prList.add(0.995);
  }


  public FindBestPrHelper(boolean debug_print){
    DEBUG_PRINT = debug_print;
  }



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

    return new double[]{bestSimuIter,bestSimuPr};
  }
  private double[] evaluatePr(int maxMemoryNum,double Pr,int succN,int failN,int[] succComNum,int[]failComNum){
    double[] simulateResult= new double[3];
    double[] successResult = simulateIteration(Pr,Pr,1,succN,maxMemoryNum,succComNum);
    double[] failResult = simulateIteration(0*(1-Pr),Pr,1,failN,maxMemoryNum,failComNum);
    simulateResult[0] = Pr*successResult[0]+(1-Pr)*(1+failResult[0]);
    simulateResult[1] = successResult[1];
    simulateResult[2] = failResult[1];
//        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\t\t\tcntPR:"+Pr+"\tsuccessN:\t"+succN+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateResult[0]+"\tsimu_nextSuccessPr:"+simulateResult[1]);
    return simulateResult;
  }
  public double findBestPr(KLLSketchLazyExact sketch, int maxMemoryNum,long rk1,long rk2,double[] deterministic_result){
    double bestEstiNum=1e9,nextSuccessPr=1.0;
    int bestPrId = 0;
    int[] successN = new int[prList.size()],failN = new int[prList.size()];
    int[][] succComNum= new int[prList.size()][],failComNum=new int[prList.size()][];
    for(int i=0;i<prList.size();i++){
      double pr=prList.getDouble(i);
      double[] cntResult = sketch.findResultRange(rk1,rk2,pr);
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
      double[] cntPrResult = evaluatePr(maxMemoryNum, prList.getDouble(i), successN[i],failN[i],succComNum[i],failComNum[i]);
      if(cntPrResult[0]<=bestEstiNum){
        bestEstiNum = cntPrResult[0];
        bestPrId = i;
        nextSuccessPr = cntPrResult[1];
      }
//            System.out.println("\t\t\t\t\tcntPR:"+prList.getDouble(i)+"\tsuccessN:\t"+successN[i]+"\t\tfailN:\t"+failN[i]+"\t\testiIter:\t"+bestEstiNum);
    }
    if(DEBUG_PRINT)System.out.println("bestPr:"+prList.getDouble(bestPrId)+"\t\testiIter:"+bestEstiNum+"\tnextSuccessN:\t"+successN[bestPrId]+"\t\t\tnextSuccessPr:\t"+nextSuccessPr);
    return prList.getDouble(bestPrId);
  }


}
