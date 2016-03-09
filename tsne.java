import java.io.*;
import java.util.Scanner;
import java.lang.*;

import no.uib.cipr.matrix.*;

/**
 * Created by 兴杲 on 2016/2/13.
 */
class q_ret{
    double[] qdata_dv;
    double[] qdata;
}

public class tsne {
    public node[] Rand_Y(int num){
        node[] y = new node[num];
        for(int i=0;i<num;i++){
            y[i] =  new node(0.0001 * Math.random(), 0.0001 * Math.random());           //it is different to set rand y value in matlab and paper(2008)
            y[i].n = i;
        }
        return y;
    }
    public double[][] readin(String fileName, int num)throws IOException{        //readin from file
        File file = new File(fileName);
        BufferedReader reader;
        String line;
        double[][] val = new double[num][];
        try{
            System.out.println("File is being read in.");
            reader = new BufferedReader(new FileReader(file));
            int linenum=0;
            int prop_num=0;

            while((((line = reader.readLine())) != null) & (linenum < num)){
                if(linenum%10==0)
                System.out.println(linenum+"lines has been read in.");
                line = line.replaceAll("\\s+","|");
                String[] arr = line.split("\\|");
                val[linenum] = new double[arr.length];
                prop_num = arr.length;
                for(int j=0,k=0;j<arr.length;j++,k++){
                    if(arr[j].length()<=0){
                        k-=1;
                        continue;
                    }
                    val[linenum][k] = Double.valueOf(arr[j]);
                }
                linenum++;
            }
            val = this.normalize(val,num,prop_num);
            System.out.println("Data Normalized.");
        }catch (IOException e){
            System.err.println("Error found when read in."+ e);
            e.printStackTrace();
        }
        return val;
    }

    public node[] run(double[][] val, int num) throws NotConvergedException {
        double iter = 1000;

        DenseMatrix X = new DenseMatrix(val);
        int col = val[0].length;
        DenseMatrix p = this.pca(X,num,col);
        DenseMatrix d_highdm = get_distance(p);
        DenseMatrix p_highdm = d2p(d_highdm);
        DenseMatrix p_s;
        p_s = (DenseMatrix) p_highdm.add(p.transpose());
        p_s = (DenseMatrix) p_highdm.scale(0.5);
        double[] dt = p_s.getData();
        double sum = 0;
        for(Double anEle : dt) sum += anEle;
        p_s = (DenseMatrix) p.scale(1/sum);
        p = p_s;
        p_s = (DenseMatrix) p.scale(4);
        dt = p_s.getData();
        double const_KL = 0;
        for(Double anEle : dt) {
            const_KL += anEle * Math.log(anEle);
        }
        node[] y = Rand_Y(num);
        double[][] yinc = new double[num][2];
        double[][] gain = new double[num][2];
        for(int i = 0; i< num; i++){
            yinc[i][0] = 0;gain[i][0] = 1;
            yinc[i][1] = 0;gain[i][1] = 1;
        }
        double momentum = 0.5;
        for (int i = 0; i < iter; i++) {
            //compute q value
            double[][] bound = {{0,0},{1,1}};
            cube q_cube = new cube(bound);
            BH_tree q_tree = new BH_tree(num, q_cube);
            q_tree.CreatLeaf(y);
            DenseMatrix q;
            q = q_tree.add_force(y);
            double[] q_val = q.getData();
            double[] q_i = new double[q.numColumns()];
            double[][] q0 = new double[q.numColumns()][q.numColumns()];
            double[][] q_dv = new double[q.numColumns()][q.numColumns()];
            for(int j = 0; j<q.numColumns(); j++){
                if(j == 0) {
                    System.arraycopy(q_val, 1, q_i, 0, q.numColumns());
                }
                else{
                    System.arraycopy(q_val, j-1, q_i, 0, j);
                    System.arraycopy(q_val, j+1, q_i, j, q.numColumns() - j);
                }
                q_ret q_ans = q_cal(q_i);
                q_i = q_ans.qdata_dv;
                if (j == 0) {
                    q0[j][0] = 0;
                    q_dv[j][0] = 0;
                    System.arraycopy(q_i, 0, q0[j], 1, q.numColumns());
                    System.arraycopy(q_ans.qdata, 0, q_dv[j], 1, q.numColumns());
                } else {
                    System.arraycopy(q_i, 0, q0[j], j - 1, j);
                    System.arraycopy(q_i, j, q0[j], j + 1, q.numColumns() - j);
                    System.arraycopy(q_ans.qdata, 0, q_dv[j], j - 1, j);
                    System.arraycopy(q_ans.qdata, j, q_dv[j], j + 1, q.numColumns() - j);
                }
            }
            DenseMatrix q_s = new DenseMatrix(q0);
            DenseMatrix Mat_minus = (DenseMatrix) p_s.add(q_s.scale(-1));

            Double[][] ygrad = new Double[num][2];
            double y_s = 0, y_s1 = 0;
            for(int j = 0; j < num; j++){
                double[] sum_j = new double[2];
                sum_j[0] = 0; sum_j[1] = 0;
                for(int k = 0; k < num; k++){
                    if(k == j) continue;
                    sum_j[0] += Mat_minus.get(j, k) * (y[j].coor[0] - y[k].coor[0]) * q_dv[j][k];
                    sum_j[1] += Mat_minus.get(j, k) * (y[j].coor[1] - y[k].coor[1]) * q_dv[j][k];
                }
                ygrad[j][0] = 4 * sum_j[0]; ygrad[j][1] = 4 * sum_j[1];
                gain[j][0] = (Math.signum(ygrad[j][0]) == Math.signum(yinc[j][0])) ? (gain[j][1] * 0.8) : (gain[j][1] + 0.2);
                gain[j][1] = (Math.signum(ygrad[j][1]) == Math.signum(yinc[j][1])) ? (gain[j][1] * 0.8) : (gain[j][1] + 0.2);
                gain[j][0] = (gain[j][0] < 0.01) ? 0.01 : gain[j][0];
                gain[j][1] = (gain[j][1] < 0.01) ? 0.01 : gain[j][1];
                yinc[j][0] = momentum * yinc[j][0] - 500 * gain[j][0] * ygrad[j][0];
                yinc[j][1] = momentum * yinc[j][1] - 500 * gain[j][1] * ygrad[j][1];
                y[j].coor[0] += yinc[j][0]; y[j].coor[1] += yinc[j][1];
                y_s += y[j].coor[0]; y_s1 += y[j].coor[1];
            }
            y_s /= num; y_s1 /= num;
            for(int j = 0; j < num; j++){
                y[j].coor[0] -= y_s;
                y[j].coor[1] -= y_s1;
            }
            if(i > 250) momentum = 0.8;
            if(i > 100) p_s = p;
        }

        return y;
    }

    public double[][] normalize(double[][] val,int num,int col){
        double[] sum = new double[col];
        double val_max=0,val_min=0;
        double[][] ans = new double[num][col];
        for(int i = 0; i< col; i++) {
            sum[i] = 0;
            for (int j = 0; j < num; j++) {
                sum[i] += val[j][i];
                if (val[j][i] > val_max) {
                    val_max = val[j][i];
                }
                if (val[j][i] < val_min) {
                    val_min = val[j][i];
                }
            }
            sum[i] /= num;
        }
        for(int i=0;i<num;i++){
            for(int j = 0;j < col; j++){
                ans[i][j] = (val[i][j]-sum[j])/(val_max-val_min);
            }
        }
        return ans;
    }

    public DenseMatrix pca(DenseMatrix val, int row, int col) throws NotConvergedException {
        System.out.println("PCA...");
        int Csize = row>col?col:row;
        DenseMatrix C = new DenseMatrix(Csize,Csize);
        if(row>=col){
            val.transBmult(val,C);
            C = (DenseMatrix) C.scale(1/col);
        }
        else{
            val.transAmult(val,C);
        }
        EVD lambda = new EVD(Csize);
        EVD l = lambda.factor(C);
        double[] b = new double[Csize];
        b = l.getRealEigenvalues();
        DenseMatrix vec = l.getRightEigenvectors();

        double[] sum = new double[Csize];
        for(int i=0;i<Csize;i++) {
            if(i==0){
                sum[i]=b[0];
            }
            else {
                sum[i] = sum[i - 1] + b[i];
            }
        }
        int num_preserve=0;
        for(int i=0;i<Csize;i++){
            if(sum[i]/sum[Csize-1]>=0.85){
                num_preserve = i;
                break;
            }
        }
        DenseMatrix pr = new DenseMatrix(Csize,num_preserve);
        int s = (Csize<num_preserve)?Csize:num_preserve;
        for(int i = 0; i < s; i++ ){
            pr.set(i,i,1);
        }
        vec.mult(pr,pr);
        double[] lambda_preserve = new double[num_preserve];
        System.arraycopy(b,0,lambda_preserve,0,num_preserve);
        if(row<col){                                                    //didn't check index carefully
            for(int i=0;i<num_preserve;i++){
                lambda_preserve[i]=1/Math.sqrt(row*lambda_preserve[i]);
            }
            val.transAmult(pr,pr);
            DenseMatrix ft = new DenseMatrix(num_preserve,num_preserve);
            for(int i=0;i<num_preserve;i++){
                ft.set(i,i,lambda_preserve[i]);
            }
            pr.mult(ft,pr);
        }
        double[] one0 = new double[row];
        for(int i = 0;i<row;i++){
            one0[i] = 1;
        }
        DenseVector one = new DenseVector(one0);
        val.mult(one,one);
        one = one.scale(row);

        for(int i=0;i<row;i++){
            for(int j=0;j<col;j++){
                val.set(i,j,val.get(i,j)-one.get(j));
            }
        }
        DenseMatrix ans = new DenseMatrix(row,num_preserve);
        val.mult(pr,ans);

        return ans;
    }

    public DenseMatrix get_distance(DenseMatrix val){
        DenseMatrix ans = new DenseMatrix(val.numRows(),val.numRows());
        int num = val.numRows();
        double sum=0;
        for(int i=0;i<num;i++){
            for(int j=0;j<num;j++){
                sum=0;
                for(int k=0;k<val.numColumns();k++) {
                    double sub = val.get(i, k)-val.get(j,k);
                    sum += sub*sub;
                }
                ans.set(i,j,sum);
                ans.set(j,i,sum);
            }
        }
        return ans;
    }

    public DenseMatrix d2p(DenseMatrix val){
        double[] data = val.getData();
        double[][] p0 = new double[val.numRows()][val.numColumns()];
        System.out.println("Converting distance to probability.");
        double[] beta = new double[val.numColumns()];
        for (int i = 0; i < beta.length; i++) beta[i] = 1;
        for (int i = 0; i < val.numColumns(); i++) {
            double s = i*10%val.numColumns();
            if(s == 0)System.out.println((i/val.numColumns()) + "percent completed");

            double[] d_i = new double[val.numColumns()-1];
            if(i == 0) {
                System.arraycopy(data, 1, d_i, 0, val.numColumns());
            }
            else{
                System.arraycopy(data, i-1, d_i, 0, i);
                System.arraycopy(data, i+1, d_i, i, val.numColumns() - i);
            }
            H_p h = Hbeta(d_i,beta[i]);

            double Hdiff = h.H - 15;            //set to 15
            int tries = 0;
            double betamin = -0xffffffff;
            double betamax = 0xffffffff;
            while(Math.abs(Hdiff) < 1e-10 && tries < 50){               //directly copied from matlab, maybe i can find a better one
                tries++;
                if(Hdiff > 0) {
                    betamin = beta[i];
                    if (betamax == 0xffffffff) beta[i] = beta[i] * 2;
                    else beta[i] = (beta[i] + betamax) / 2;
                }
                else{
                    betamax = beta[i];
                    if(betamin == -0xffffffff) beta[i] = beta[i] / 2;
                    else beta[i] = (beta[i] + betamin) / 2;
                }
                h = Hbeta(d_i,beta[i]);
                Hdiff = h.H - 15;
            }
            if(i == 0) {
                p0[i][0] = 0;
                System.arraycopy(h.p, 0, p0[i], 1, val.numColumns());
            }
            else{
                System.arraycopy(h.p, 0, p0[i], i-1, i);
                System.arraycopy(h.p, i, p0[i], i+1, val.numColumns() - i);
            }
        }
        return new DenseMatrix(p0);
    }

    private q_ret q_cal(double[] q){
        q_ret ans = new q_ret();
        double sum = 0;
        ans.qdata = new double[q.length];
        ans.qdata_dv = new double[q.length];
        for(int i = 0; i < q.length; i++){
            ans.qdata[i] = 1/(1+q[i]);
            sum += ans.qdata[i];
        }
        for(int i = 0; i< q.length; i++){
            ans.qdata_dv[i] = ans.qdata[i] / sum;
        }
        return ans;
    }

    private H_p Hbeta(double[] d, double beta){        //slightly changed from matlab_tsne
        H_p ans = new H_p();
        ans.p = new double[d.length];
        double sumpp = 0;
        double sump = 0;
        for(int i=0;i<d.length;i++){
            ans.p[i] = Math.exp(-d[i] * beta);
            sump += ans.p[i];
            sumpp += ans.p[i]*Math.log(ans.p[i]);
        }
        ans.H = Math.log(sump) - sumpp/sump;
        for (int i = 0; i < d.length; i++) {
            ans.p[i] = ans.p[i] / sump;
        }
        return  ans;
    }

    public static void main(String[] arg) throws IOException, NotConvergedException {
        tsne Tsne = new tsne();
        double[][] val;
        Scanner rd = new Scanner(System.in);
        System.out.println("Filename:");
        String filename = rd.nextLine();
        System.out.println("Number of samples:");
        int Sample_num = rd.nextInt();
        val = Tsne.readin(filename,Sample_num);
        node[] ans = Tsne.run(val, Sample_num);
        System.out.println("Filename to write in:");
        String wt = rd.nextLine();
        try{
            FileWriter fileWriter = new FileWriter(wt);
            String ln;
            for(int i = 0; i < Sample_num; i++) {
                ln = ans[i].coor[0] + " " + ans[i].coor[1] + "\n";
                fileWriter.write(ln);
            }
            fileWriter.close();
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }
}

class H_p{
    double H;
    double[] p;
    H_p(){H=0;}
}
