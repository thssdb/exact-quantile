import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.util.Random;

public class MainForCombinationNumber {

    public static void main(String[] args) {
//        NormalDistribution dis = new NormalDistribution(0, 1.0);
//        System.out.println("\t\t\ts:"+dis.cumulativeProbability(-0.5));
//        BinomialDistribution dis2 = new BinomialDistribution(100/2,0.5);

        double[][] pr = new double[5000][5000];
        pr[1][0]=pr[1][2]=0.25;pr[1][1]=0.5;
        for(int m=2;m<=2333;m++){
            double tmpSum=0;
            for(int sum=-m;sum<=m;sum++) {
                if(m - 1 + sum - 1>=0)pr[m][m + sum] += 0.25 *pr[m - 1][m - 1 + sum - 1];
                if(m - 1 + sum>=0)pr[m][m + sum] += 0.5 * pr[m - 1][m - 1 + sum];
                pr[m][m + sum] += 0.25 * pr[m - 1][m - 1 + sum + 1];

                tmpSum+=pr[m][m + sum];
            }
//            System.out.println("\t\t\tsum="+tmpSum);
//            for(int sum=-m;sum<=m;sum++)System.out.print("\t"+pr[m][m+sum]*Math.pow(2.0,m*2));System.out.println();
//            for(int sum=-m;sum<=m;sum++)System.out.print("\t"+pr[m][m+sum]);System.out.println();
        }
        int m=32,M=30;
        for(int sum=-M;sum<=M;sum++)System.out.print("\t"+pr[m][m+sum]);System.out.println();
        double mu = 0,sig = Math.sqrt(m/2.0);
        DoubleArrayList tmp = new DoubleArrayList();
        Random random = new Random();
        int RANDOM_CASE=10000000;
        for(int i=0;i<RANDOM_CASE;i++)tmp.add((random.nextGaussian())*sig);
        double mxDeltaRate=0;
        for(int i=-M;i<=M;i++)if(pr[m][m+i]>=1e-3){System.out.println("\t\t\t"+i);
            int cnt=0;
            for(double j:tmp)
                if(j>=i-0.5&&j<=i+0.5)cnt++;
//            System.out.print("\t"+1.0*cnt/RANDOM_CASE);
//            System.out.print("\t"+(1.0*cnt/RANDOM_CASE)/pr[m][m+i]);
            double delRate = (1.0*cnt/RANDOM_CASE)/pr[m][m+i];
            mxDeltaRate = Math.max(mxDeltaRate,delRate);
        }System.out.println("\n"+"\tMxDelRate:"+mxDeltaRate);
    }
}
