/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.geometry.contourtree

/**
 *
 * Join tree has its root at the top
 * Split tree has its root at the bottom
 */
class JSTreeNode(val vertex: Int) {
  var joinParent: JSTreeNode = null
  var numJoinChildren: Int = 0

  var splitParent: JSTreeNode = null
  var numSplitChildren: Int = 0
  var deleted: Boolean = false
}