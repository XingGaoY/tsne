/**
 * Created by 兴杲 on 2016/2/19.
 */
import java.util.*;
import no.uib.cipr.matrix.DenseMatrix;

class new_bound{
    double[][] bound;
    int ind;
    new_bound(double[][] _bound,int _ind){
        ind = _ind;
        bound = _bound;
    }
}
class node{
    double[] coor;
    int mass;
    int n;
    node(){}
    node(double _x,double _y){
        coor = new double[2];
        coor[0] = _x;
        coor[1] = _y;
        mass = 1;
    }
    public new_bound contain(int axis, new_bound _bound){
        double mid = (_bound.bound[axis-1][1]+_bound.bound[axis-1][0])/2;
        if(coor[axis-1]>mid){
            _bound.ind |= axis;
            _bound.bound[axis-1][0] = mid;
        }
        else
        {
            _bound.bound[axis-1][1] = mid;
        }
        return _bound;
    }

    public double compute_force(node p){
        double x = coor[0] - p.coor[0];
        double y = coor[1] - p.coor[1];
        return x*x + y*y;
    }

}
class cube{
    double[][] location;
    cube(double[][] _loc)
    {
        location = _loc;
    }
    boolean is_in (node p){
        return(p.coor[0]>location[0][0]&&
                p.coor[0]<=location[0][1]&&
                p.coor[1]>location[1][0]&&
                p.coor[1]<=location[1][1]);
    }
}
public class BH_tree {
    cube Cube;
    node Cell;
    BH_tree[] child;
    Integer[] ele;

    BH_tree(){
        Cell.coor =new double[2];
        Cell.coor[0] = 0;
        Cell.coor[1] = 0;
        Cell.mass = 0;
    }
    BH_tree(int num,cube _bound){
        for(int i=0;i<num;i++) ele[i] = i;
        Cube = _bound;
        Cell = new node();
        Cell.mass = num;
        Cell.coor =new double[2];
        Cell.coor[0] = 0.5;
        Cell.coor[1] = 0.5;
    }

    public void CreatLeaf(node[] y){
        Vector[] new_ele = new Vector[8];
        if(ele.length == 1){
            child = null;           //leaf Node
        }
        else{
            child = new BH_tree[4];
            for (Integer anEle : ele) {
                new_bound compute_bound = new new_bound(Cube.location, 0);
                compute_bound = y[anEle].contain(1, compute_bound);
                compute_bound = y[anEle].contain(2, compute_bound);
                if (child[compute_bound.ind] == null) {
                    child[compute_bound.ind] = new BH_tree();
                }
                child[compute_bound.ind].Cube = new cube(compute_bound.bound);
                child[compute_bound.ind].Cell.mass += 1;
                child[compute_bound.ind].Cell.coor[0] += y[anEle].coor[0];
                child[compute_bound.ind].Cell.coor[1] += y[anEle].coor[1];
                new_ele[compute_bound.ind].add(anEle);
            }
            for(int i=0;i<4;i++){
                child[i].ele = (Integer[]) new_ele[i].toArray();
                child[i].Cell.coor[0] /= child[i].ele.length;
                child[i].Cell.coor[1] /= child[i].ele.length;
                child[i].CreatLeaf(y);
            }
        }
    }

    public DenseMatrix add_force(node[] y){                          //compute a distance twice is a waste.
        DenseMatrix ans = new DenseMatrix(y.length,y.length);
        LinkedList<BH_tree> list = new LinkedList<BH_tree>();
        for(int i=0; i<y.length; i++) {
            node p = y[i];
            BH_tree b = this;
            if (b != null) list.add(b);
            while (list.isEmpty()) {
                b = list.getFirst();
                list.removeFirst();
                if (!b.Cube.is_in(p)) {
                    double d = Math.sqrt(b.Cell.compute_force(p));
                    double s = Math.max(b.Cube.location[0][1] - b.Cube.location[0][0], b.Cube.location[1][1] - b.Cube.location[1][0]);
                    double theta = 0.5;
                    if (s / d < theta && b.ele.length != 1) {
                        for (Integer anEle : ele) {
                            ans.set(p.n, b.ele[anEle], d * d);
                            ans.set(b.ele[anEle], p.n, d * d);
                        }
                    } else {
                        for (int j = 0; j < 4; j++) {
                            if (b.child[j] != null)
                                list.addLast(b.child[j]);
                        }
                    }
                }
            }
        }
        return ans;
    }
}
