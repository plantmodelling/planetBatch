/* Created on May 30, 2005 */
package maize;

import javax.vecmath.Tuple3f;

/** @author Xavier Draye - Universit� catholique de Louvain (Belgium) */
public class SpaceTimeInfo {
   public double x, y, z, t;
   public String code;
   private static long nextFreeCode = 1L; 
   
   public SpaceTimeInfo(String code) {
      this.code = code;
   }
   
   public SpaceTimeInfo(float x, float y, float z, float t, String code) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.t = t;
      this.code = code;
   }
   
   public long getNewCode() {
      return nextFreeCode++;
   }
   
   public SpaceTimeInfo setCoordinates(Tuple3f c, float t) {
      x = c.x;
      y = c.y;
      z = c.z;
      this.t = t;
      return this;
   }

   public SpaceTimeInfo setCoordinates(float x, float y, float z, float t) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.t = t;
      return this;
   }

   public SpaceTimeInfo setCoordinates(Tuple3f c) {
      x = c.x;
      y = c.y;
      z = c.z;
      return this;
   }
      
   public SpaceTimeInfo setCoordinates(float x, float y, float z) {
      this.x = x;
      this.y = y;
      this.z = z;
      return this;
   }
   
   public SpaceTimeInfo setCoordinates(double x, double y, double z) {
	      this.x = x;
	      this.y = y;
	      this.z = z;
	      return this;
	   }
      
   public SpaceTimeInfo setX(float x) {
      this.x = x;
      return this;
   }
   
   public SpaceTimeInfo setY(float y) {
      this.y = y;
      return this;
   }
   
   public SpaceTimeInfo setZ(float z) {
      this.z = z;
      return this;
   }
   
   public SpaceTimeInfo setT(float t) {
      this.t = t;
      return this;
   }
   
   
   public float getX() {return (float) x;}
   public double getXd() {return x;}
   public float getY() {return (float) y;}
   public double getYd() {return x;}
   public float getZ() {return (float) z;}
   public double getZd() {return x;}
   public float getT() {return (float) t;}
   public double getTd() {return x;}
}
