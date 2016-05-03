/*
Class to generate a .dat file for faces and lines.

The generated file can be loaded using ZBuf.java or Hlines.java

However, The ZBuf.java won't show the climbable supports.

The program can be run using command Line as follows

java Stairs 25 30 stairs.dat

where 25 is the number of stairs; 30 is the angle of separation between two stairs and stairs.dat is the file where the vertices and faces will be defined.
 */


/**
 *
 * @author Dhananjay Singh
 * @Net Id: dxs145530
 * @UTD Id: 2021250625
 */
import java.io.*;
public class Stairs 
{
   int n; //This is the number of faces of the central Cylinder
   int sn; //No. of stairs.
   double alphaDeg; //angle of seperation between stairs
   String fName; //Output file
   FileWriter fw;
   float rOuter; //inner and outer radius of the central cylinder
   float rInner;
   public Stairs(int sn,double deg,String fName) throws IOException
   {  
      n=50;
      this.sn=sn;
      alphaDeg=deg;
      this.fName=fName;
      fw = new FileWriter(fName);
      rOuter=2.5f/2.0f;
      rInner=0.0f;
   }
   void genStair()throws IOException
   {  
      // To generate stairs .
      int t=0;
      double a = 7;
      
      //The following code is reused directly from Beams.java
      alphaDeg=alphaDeg * Math.PI / 180;
      Point3D[] P = new Point3D[11];
      P[1] = new Point3D(a, -1, 0);
      P[2] = new Point3D(a,  1, 0);
      P[3] = new Point3D(1,  1, 0);
      P[4] = new Point3D(1, -1, 0);
      P[5] = new Point3D(a, -1, 1);
      P[6] = new Point3D(a,  1, 1);
      P[7] = new Point3D(1,  1, 1);
      P[8] = new Point3D(1, -1, 1);
      
      //The support beams
      P[9]=new Point3D(7.0,0.0,0.1);
      P[10]=new Point3D(7.0,0.0,6.0);
      for (int k=0; k<sn; k++)
      {  // Beam k:
         double phi = k * alphaDeg,
            cosPhi = Math.cos(phi), sinPhi = Math.sin(phi);
         int m = 10 * k;
         for (int i=1; i<=10; i++)
         {  double x = P[i].x, y = P[i].y;
            float x1 = (float)(x * cosPhi - y * sinPhi),
                  y1 = (float)(x * sinPhi + y * cosPhi),
                  z1 = (float)(P[i].z + k);
            fw.write(( m + i) + " " + x1 + " " + y1 + " " + z1 + 
               "\r\n");
            t=(m+i)+1; //Keeping track of vertex number as Cylinder.java and Beam.java are working in unison here
         }
      }
      //Generate vertices and faces of the cylinder ,dependent on the vertices of the last generated beam
      generateCylinder(t,rOuter,rInner,n);
      
      //Generating faces for the stairs and support beams
      for (int k2=0; k2<sn; k2++)
      {  // Beam k again:
         int m = 10 * k2;
         face(m, 1, 2, 6, 5);
         face(m, 4, 8, 7, 3);
         face(m, 5, 6, 7, 8);
         face(m, 1, 4, 3, 2);
         face(m, 2, 3, 7, 6);
         face(m, 1, 5, 8, 4);
         //To generate the support beam faces (a line in this case
         face1(m,9,10,sn);
         //Note that the two beams should be connected to one another
         face1(m,10,20,sn);
      }
      fw.close();
}
void face(int m, int a, int b, int c, int d)throws IOException
   {  a += m; b += m; c += m; d += m;
      fw.write(a + " " + b + " " + c + " " + d + ".\r\n");
   }
void face1(int m,int a,int b,int sn) throws IOException
{
    a+=m;
    b+=m;
    
    //The last support horizontal beam should not connect to the last step, so we perform the following check
    if(b<=10*sn)
    fw.write(a+" "+b+".\r\n");
}

//Re used directly from Cylinder.java
void generateCylinder(int vertex,float rOuter,float rInner ,int faces) throws IOException
{
    int n2 = 2 * faces, n3 = 3 * faces, n4 = 4 * faces;
    double delta = 2 * Math.PI / faces;
    for (int i=1; i<=faces; i++)
      {  double alpha = i * delta,
            cosa = Math.cos(alpha), sina = Math.sin(alpha);
         for (int inner=0; inner<2; inner++) 
         {  double r = (inner == 0 ? rOuter : rInner);
            if (r > 0) 
            for (int bottom=0; bottom<2; bottom++) 
            {  int k = (2 * inner + bottom) * faces + i+vertex;
               // Vertex numbers for i = 1:
               // Top:       1 (outer)   2n+1 (inner)
               // Bottom:  n+1 (outer)   3n+1 (inner)
               wi(k); // w = write, i = int, r = real
               wr(r * cosa); wr(r * sina); // x and y
               wi((faces/2)-((faces/2)*bottom)); // bottom: z = 0; top: z = 25
               fw.write("\r\n");
            }
         }
      }
    fw.write("Faces:\r\n");
    // Top boundary face:
      for (int i=1; i<=faces; i++) wi(i+vertex);
      if (rInner > 0)
      {  wi(-1*(n3+vertex)); // Invisible edge, see Section 7.5
         for (int i=n3-1; i>=n2+1; i--) wi(i+vertex);
         wi(n3+vertex); wi(-1*(faces+vertex)); // Invisible edge again.
      }
      fw.write(".\r\n");
      // Bottom boundary face:
      for (int i=n2; i>=faces+1; i--) wi(i+vertex);
      if (rInner > 0)
      {  wi(-1*(n3+1+vertex));
         for (int i=n3+2; i<=n4; i++) wi(i+vertex);
         wi(n3+1+vertex); wi(-1*(faces+1+vertex));
      }
      fw.write(".\r\n");
      // Vertical, rectangular faces:
      for (int i=1; i<=faces; i++)
      {  int j = i % faces + 1;
         // Outer rectangle:
         wi(j+vertex); wi(i+vertex); wi(i + faces+vertex); wi(j + vertex+faces); fw.write(".\r\n");
         if (rInner > 0)
         {  // Inner rectangle:
            wi(i + n2 + vertex); wi(j + n2+vertex); wi(j + n3+vertex); wi(i + n3+vertex); 
            fw.write(".\r\n"); 
         }
      }
}
   void wi(int x) throws IOException
   {  fw.write(" "+String.valueOf(x));
   }
   
   void wr(double x) throws IOException
   {  if (Math.abs(x) < 1e-9) x = 0;
      fw.write(" "+String.valueOf((float)x));
      // float instead of double to reduce the file size
   }
   
   //Entry point of the program
   // If the user does not pass any command line arguments, default values will be used
   // The default output file is TestStairs.dat unless specified by the user
   public static void main(String args[])
   {
       int n=25;
       double degrees=30.0d;
       String fName="TestStairs.dat";
       for(String s:args)
       {
           n=Integer.parseInt(args[0]);
           degrees=Double.parseDouble(args[1]);
           fName=args[2];
       }
       try
       {
       Stairs stairs=new Stairs(n,degrees,fName);
       stairs.genStair();
       }catch(Exception e){}
   }
}
