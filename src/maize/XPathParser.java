/* Created on Feb 18, 2006 */
package maize;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;
//import javax.xml.xpath.XPathFactoryConfigurationException;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;


/** XPathParser provides direct access to the content of an XML document by means of
 * XPath expressions. An XPathParser object is initialized with an expression that
 * sets a context node and allows the extraction of String, int, float, boolean
 * (single items, or arrays, or array items) at locations relative to the context node.</p>
 * The class provides constructors that take the filename of the XML file and the expression
 * of the context node, as well as constructors that take another XPathParser and a (relative)
 * expression of the context node. This architecture allows easy traversal of the XML tree.</p>
 * 
 * Let us suppose the following part of an XML file:</p>
 * 
 * <cnode>
 *  <var0>Hallo</var0>
 *  <var1>1.0</var1>
 *  <var2>true</var2>
 *  <vector3>1</vector3>
 *  <vector3>2</vector3>
 *  <vector3>3</vector3>
 *  <vector3>4</vector3>
 *  <matrixrow4>
 *   <vector4>m1.1</vector4>
 *   <vector4>m1.2</vector4>
 *  </var4>
 *  <matrixrow4>
 *   <vector4>m2.1</vector4>
 *   <vector4>m2.2</vector4>
 *  </var4>
 * </cnode>
 * 
 * </p>
 * The following code retreives all data from the XML file:</p>
 * <code>
 * // Open the XML file and set the expression that points to cnode (the context node)</p>
 * XPathParser xpp = new XPathParser("file.xml", "/cnode");</p>
 * </p>
 * // Retreive the values of elements just below the context node</p>
 * String var0 = xpp.getString("Hallo");</p>
 * float var1 = xpp.getFloat("var1");</p>
 * boolean var2 = xpp.getBoolean("var2");</p>
 * </p>
 * // Find out the number of occurences of element vector3</p>
 * int n = xpp.getElementCount("vector3");</p>
 *  </p>
 * // Elements that appear more than once (ex: a vector) can be accessed one by one</p>
 * int vector31 = xpp.getIntegerItem("vector3", 1);</p>
 * int vector32 = xpp.getIntegerItem("vector3", 2);</p>
 * int vector33 = xpp.getIntegerItem("vector3", 3);</p>
 * int vector34 = xpp.getIntegerItem("vector3", 4);</p>
 * </p>
 * // Elements that appear more than once can also be returned at once using arrays</p>
 * int[] vector3 = xpp.getIntegerArray("vector3");</p>
 * </p>
 * // Collections of elements that appear more than once can be retreived by creating</p>
 * // new XPathParsers for each collection and getting elements under each. We use a </p>
 * // specific constructor that creates an XPathParser whose context node is defined</p>
 * // relative to the current XPathParser and considering its rank in the set of</p>
 * // collections. Example: a matrix has a number of row vectors, each being a collection</p>
 * // of elements.</p>
 * int n = xpp.getElementCount("matrixrow4");</p>
 * Vector matrix4 = new Vector();</p>
 * for (int i = 0; i < n; i++) {</p>
 *    XPathParser xpp1 = new XPathParser(xpp, "matrix4", i);</p>
 *    String[] vector4 = xpp1.getStringArray("vector4");</p>
 *    matrix4.add(vector4);</p>
 * }</p>
 * </code>
 *  @author Xavier Draye - Universit� catholique de Louvain (Belgium)
 */
public class XPathParser {

   String contextNode;
   // Either an InputSource or an Element
   Object source;
   static XPath xp = XPathFactory.newInstance().newXPath();
   static final String defaultContextNode = ".";
   
   // TODO Can discard this - not necessary since 1.5.0_06
   // Get a unique instance of XPath that will be shared by all XPathParser instances
//   {
//      try {
//         xp = XPathFactory.newInstance().newXPath();
//      } catch (XPathFactoryConfigurationException e) {
//         Util.log("XML error: can't get an XPath.");
//         Util.logException(e);
//      }
//   }
   
   /** Creates a new instance of XPathParser from an InputSource or Element and the contextNode
    *  @param o An InputSource pointing to the XML file or an Element in an existing Document instance 
    *  @param contextNode An XPath expression which sets the absolute location of the 
    * context node (therefore starting with "/")
    */ 
   public XPathParser(Object o, String contextNode) {
      this.source = o;
      this.contextNode = contextNode;
   }
   
   /** Creates a new instance of XPathParser from an InputSource or Element. The context node is
    * set to the root of the InputSource or to the Element.
    *  @param o An InputSource pointing to the XML file or an Element in an existing Document instance 
    * context node (therefore starting with "/")
    */ 
   public XPathParser(Object o) {
      this.source = o;
      this.contextNode = defaultContextNode;
   }
   
   /** Creates a new instance of XPathParser from the XML filename and the contextNode
    *  @param xmlFileName The filename of the XML file 
    *  @param contextNode An XPath expression which sets the absolute location of the 
    * context node (therefore starting with "/")
    */ 
   public XPathParser(String xmlFileName, String contextNode) {
      this(new InputSource(xmlFileName), contextNode);
   }
   
   /** Creates a new instance of XPathParser from the XML filename and set the contextNode 
    * to the root of the document ("/")
    *  @param xmlFileName The filename of the XML file 
    */ 
   public XPathParser(String xmlFileName) {
      this(new InputSource(xmlFileName), defaultContextNode);
   }


   /** Creates a new instance of XPathParser relative to an existing XPathParser object. If
    * the node set of the XPathExpression contains more than one element, the context node of
    * the new object is the first of those elements
    *  @param xpr The existing XPathParser relative to which the context node is defined 
    *  @param subExpression An XPath expression which sets the location of the 
    * context node relative to the context node of the existing XPathParser.
    */ 
   public XPathParser(XPathParser xpr, String subExpression) {
      this(xpr.getSource(), xpr.getContextNode() + "/" + subExpression);
   }

   /** Creates a new instance of XPathParser relative to an existing XPathParser object.
    * Use this constructor when the XPathExpression contains more than one element and the
    * context node of the new object should not be the first element 
    *  @param xpr The existing XPathParser relative to which the context node is defined 
    *  @param subExpression An XPath expression which sets the location of the 
    * context node relative to the context node of the existing XPathParser
    *  @param i The rank of the element that should serve as context node
    */ 
   public XPathParser(XPathParser xpr, String subExpression, int i) {
      this(xpr.getSource(), xpr.getContextNode() + "/" + subExpression + "[" + (i + 1) + "]");
   }

   public String getContextNode() {return contextNode;}

   public void setContextNode(String s) {
      contextNode = contextNode.equals("") ? defaultContextNode : s;
      }
   
   public Object getSource() {return source;}
   
   public void setSource(Object o) {this.source = o;}
   
   /** Retreive a double value from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the number will be retreived from.
    *  @return A double value.
    */
   public double getNumber(String subExpression) {
      double result = 0.0;
      String expression = contextNode + "/" + subExpression;
      try {
         result = (Double)xp.evaluate(expression, source, XPathConstants.NUMBER);
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }
   
   /** Retreive an int value from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the int will be retreived from.
    *  @return An int value.
    */
   public int getInteger(String subExpression) {
      return (int)getNumber(subExpression);
   }

   /** Retreive an int value from its location relative to the context node. To be used 
    * when the XPathExpression results in more than one elements and the int value should
    * not be retreived from the first element 
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the int will be retreived from
    *  @param i The rank of the element where the int must be retreived from 
    *  @return An int value.
    */
   public int getIntegerItem(String subExpression, int i) {
      int result = 0;
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         if (i < nl.getLength()) result = Integer.parseInt(nl.item(i).getTextContent());
      } catch(Exception e) {
         Util.log("XML error: can't read node " + subExpression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a float value from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the float will be retreived from.
    *  @return A float value.
    */
   public float getFloat(String subExpression) {
      return (float)getNumber(subExpression);
   }
   
   /** @see #getIntegerItem(String, int)
    */
   public float getFloatItem(String subExpression, int i) {
      float result = 0f;
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         if (i < nl.getLength()) result = Float.parseFloat(nl.item(i).getTextContent());
      } catch(Exception e) {
         Util.log("XML error: can't read node " + subExpression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a String value from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the String will be retreived from.
    *  @return A String.
    */
   public String getString(String subExpression) {
      String result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         result = (String)xp.evaluate(expression, source, XPathConstants.STRING);
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.log(e.getMessage());
      }
      return result;
   }

   /** @see #getIntegerItem(String, int)
    */
   public String getStringItem(String subExpression, int i) {
      String result = "";
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         if (i < nl.getLength()) result = nl.item(i).getTextContent();
      } catch(Exception e) {
         Util.log("XML error: can't read node " + subExpression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a boolean value from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the boolean will be retreived from.
    *  @return A boolean.
    */
   public boolean getBoolean(String subExpression) {
      boolean result = false;
      String expression = contextNode + "/" + subExpression;
      try {
         result = (Boolean) xp.evaluate(expression, source, XPathConstants.BOOLEAN);
         Util.log("BOOLEAN = "+result+ " / source = "+expression);
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** @see #getIntegerItem(String, int)
    */
   public boolean getBooleanItem(String subExpression, int i) {
      boolean result = false;
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         if (i < nl.getLength()) result = Boolean.parseBoolean(nl.item(i).getTextContent());
      } catch(Exception e) {
         Util.log("XML error: can't read node " + subExpression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a ParameterSet object from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the ParameterSet object will be retreived from.
    *  @return A ParameterSet object.
    */
   public ParameterSet getParameterSet(String subExpression) {
      ParameterSet result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         Node n = (Node)xp.evaluate(expression, source, XPathConstants.NODE);
         result = (ParameterSet)n.getUserData("parameterSetObject");
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** @see #getIntegerItem(String, int)
    */
   public ParameterSet getParameterSetItem(String subExpression, int i) {
      ParameterSet result = null;
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         if (i < nl.getLength()) result = (ParameterSet)nl.item(i).getUserData("parameterSetObject");
      } catch(Exception e) {
         Util.log("XML error: can't read node " + subExpression + ".");
         Util.logException(e);
      }
      return result;
   }


   /**
    * Check whether the DOM Element corresponding to subExpression is a template element. A
    * template element, say te, is one that is present to tell the user that te elements
    * can can optionally be placed there. Template elements should be skipped when a ParameterSet
    * object is parsing the DOM to read its parameter values/ 
    * @param subExpression The XPath expression (relative to the context node) of the element
    * @return true if subExpression points to a template element
    */
   public boolean isTemplate(String subExpression) {
      Element e = getElement(subExpression);
      if (e == null) return false;
      Boolean isTemplate = (Boolean) e.getUserData("isTemplate");
      return isTemplate == null ? false : isTemplate;
   }
   

   /** Retreive an Element from its location relative to the context node
    *  @param subExpression The XPath expression (relative to the context node) where
    *  the Element will be retreived from.
    *  @return A ParameterSet object.
    */
   public Element getElement(String subExpression) {
      Element result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         Node n = (Node)xp.evaluate(expression, source, XPathConstants.NODE);
         if (n.getNodeType() == Node.ELEMENT_NODE) result = (Element) n;
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }
   

   /** @see #getIntegerItem(String, int)
    */
   public Element getElementItem(String subExpression, int i) {
      Element result = null;
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         Node n = null;
         if (i < nl.getLength()) n = nl.item(i);
         if (n != null && n.getNodeType() == Node.ELEMENT_NODE) result = (Element) n;
      } catch(Exception e) {
         Util.log("XML error: can't read node " + subExpression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Returns the size of the node set defined by an XPath expression relative to
    * the context node, or in other words, the number of occurences of the element
    * defined by the expression. 
    *  @param subExpression The XPath expression (relative to the context node) to be evaluated
    *  @return The number of occurences. The method returns zero if there are no elements, if
    *  subexpression is invalid or is there is a single element and this element is template 
    *  (@see isTemplate(String))
    */
   public int getElementCount(String subExpression) {
      int result = 0;
      String expression = contextNode + "/" + subExpression;
      try {
         NodeList nl = (NodeList)xp.evaluate(expression, source, XPathConstants.NODESET);
         result = nl.getLength();
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      if ((result == 1) && isTemplate(subExpression)) result = 0;
      return result;
   }
   
   
   /** Retreive a collection of strings from their location relative to the context node.
    *  @param subExpression The XPath expression (relative to the context node) to be evaluated
    *  @return A new string array whose size is the number of elements resulting from the XPath expression
    */
   public String[] getStringArray(String subExpression) {
      String[] result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         NodeList nl = (NodeList)xp.evaluate(contextNode + "/" + subExpression, source, XPathConstants.NODESET);
         int n = nl.getLength();
         if (n > 0) {
            result = new String[n];
            for (int i = 0; i < n; i++) result[i] = nl.item(i).getTextContent();
         }
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a collection of integer from their location relative to the context node.
    *  @param subExpression The XPath expression (relative to the context node) to be evaluated
    *  @return A new int array whose size is the number of elements resulting from the XPath expression
    */
   public int[] getIntegerArray(String subExpression) {
      int[] result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         NodeList nl = (NodeList)xp.evaluate(expression, source, XPathConstants.NODESET);
         int n = nl.getLength();
         if (n > 0) {
            result = new int[n];
            for (int i = 0; i < n; i++) result[i] = Integer.parseInt(nl.item(i).getTextContent());
         }
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a collection of float from their location relative to the context node.
    *  @param subExpression The XPath expression (relative to the context node) to be evaluated
    *  @return A new float array whose size is the number of elements resulting from the XPath expression
    */
   public float[] getFloatArray(String subExpression) {
      float[] result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         NodeList nl = (NodeList)xp.evaluate(expression, source, XPathConstants.NODESET);
         int n = nl.getLength();
         if (n > 0) {
            result = new float[n];
            for (int i = 0; i < n; i++) result[i] = Float.parseFloat(nl.item(i).getTextContent());
         }
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }

   /** Retreive a collection of boolean from their location relative to the context node.
    *  @param subExpression The XPath expression (relative to the context node) to be evaluated
    *  @return A new boolean array whose size is the number of elements resulting from the XPath expression
    */
   public boolean[] getBooleanArray(String subExpression) {
      boolean[] result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         NodeList nl = (NodeList)xp.evaluate(expression, source, XPathConstants.NODESET);
         int n = nl.getLength();
         if (n > 0) {
            result = new boolean[n];
            for (int i = 0; i < n; i++) result[i] = Boolean.parseBoolean(nl.item(i).getTextContent());
         }
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }


   /** Retreive a collection of ParameterSets from their location relative to the context node.
    *  @param subExpression The XPath expression (relative to the context node) to be evaluated
    *  @return A new ParameterSet array whose size is the number of elements resulting from the XPath expression
    */
   public ParameterSet[] getParameterSetArray(String subExpression) {
      ParameterSet[] result = null;
      String expression = contextNode + "/" + subExpression;
      try {
         NodeList nl = (NodeList)xp.evaluate(expression, source, XPathConstants.NODESET);
         int n = nl.getLength();
         // If there is a single item and it is a template, drop it 
         if (n == 1 && (Boolean)nl.item(0).getUserData("isTemplate")) 
            return new ParameterSet[0];
         if (n > 0) {
            result = new ParameterSet[n];
            for (int i = 0; i < n; i++) {
               result[i] = (ParameterSet) nl.item(i).getUserData("parameterSetObject");
            }
         }
      } catch(Exception e) {
         Util.log("XML error: can't read node " + expression + ".");
         Util.logException(e);
      }
      return result;
   }


}
