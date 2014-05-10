/**
 * Ayla Geometry Core
 * (c) 2011-2014 William Harvey
 * http://www.aylasoftware.org
 */
package ayla.util

object Timing {
	def apply[R](task: String)(f: => R): R = {
		println("Performing task: " + task)
		val t1 = System.currentTimeMillis
		val result = f
		val t2 = System.currentTimeMillis
		val seconds = (t2 - t1) / 1000.0
		println("Time(" + task + "): " + seconds + " seconds.")
		return result
	}
}
