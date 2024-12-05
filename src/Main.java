import java.util.HashMap;

public class Main {
    public static void main(String[] args) {
        Function f = new Function();
        if (args[0].equals("gwas")){
            f.fbp(args[1]);
        }else{
            System.out.println("Invaid option " + args[0]);
        }
    }
}