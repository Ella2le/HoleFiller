import java.io.File

data class Vertex(val id: Int, val x: Float, val y: Float, val z: Float) {
    val adjacentVertices = mutableListOf<Int>()
    val adjacentEdges = mutableListOf<Int>()
    val adjacentTriangles = mutableListOf<Int>()
}

data class Edge(val id: Int, val v1: Int, val v2: Int)

data class Triangle(val id: Int, val v1: Int, val v2: Int, val v3: Int)

class Mesh {
    val triangles = mutableListOf<Triangle>()
    val vertices = mutableListOf<Vertex>()
    var edges = mutableListOf<Edge>()

    fun addVertex(x: Float, y: Float, z: Float) : Int {
        val id = if (vertices.isEmpty()) 0 else vertices.last().id + 1
        vertices.add(Vertex(id, x, y, z))
        return id
    }

    fun makeVerticesNeighbors(v1: Int, v2: Int) : Boolean {
        if (vertices[v1].adjacentVertices.contains(v2))
            return false

        vertices[v1].adjacentVertices.add(v2)
        vertices[v2].adjacentVertices.add(v1)
        return true
    }

    fun makeVerticesUnneighbors(v1: Int, v2: Int) : Boolean {
        if(!vertices[v1].adjacentVertices.contains(v2))
            return false

        vertices[v1].adjacentVertices.remove(v2)
        vertices[v2].adjacentVertices.remove(v1)
        return true
    }

    fun addEdge(v1: Int, v2: Int) {
        val id= if (edges.isEmpty()) 0 else edges.last().id + 1
        edges.add(Edge(id, v1, v2))
        vertices[v1].adjacentEdges.add(id)
        vertices[v2].adjacentEdges.add(id)
    }

    fun removeEdge(v1: Int, v2: Int) {
        var id = -1
        for(edge in edges) {
            if((edge.v1 == v1 && edge.v2 == v2) || (edge.v2 == v1 && edge.v1 == v2)) {
                id = edge.id
                edges.remove(edge)
            }
        }
        if(id == -1) return
        vertices[v1].adjacentVertices.remove(id)
        vertices[v2].adjacentVertices.remove(id)
    }

    fun triangleExists(v1: Int, v2: Int, v3: Int) : Boolean {
        for(triangle in triangles) {
            if(triangle.v1 == v1 && triangle.v2 == v2 && triangle.v3 == v3) {
                return true
            }
        }
        return false
    }

    fun addTriangle(v1: Int, v2: Int, v3: Int) {
        val id = if (triangles.size == 0) 0 else triangles.last().id + 1
        triangles.add(Triangle(id, v1, v2, v3))

        vertices[v1].adjacentTriangles.add(id)
        vertices[v2].adjacentTriangles.add(id)
        vertices[v3].adjacentTriangles.add(id)

        if(makeVerticesNeighbors(v1, v2)) addEdge(v1, v2)
        if(makeVerticesNeighbors(v1, v3)) addEdge(v1, v3)
        if(makeVerticesNeighbors(v2, v3)) addEdge(v2, v3)
    }

    fun removeTriangle(v1: Int, v2: Int, v3: Int) {
        var id = -1
        for(triangle in triangles) {
            if(triangle.v1 == v1 && triangle.v2 == v2 && triangle.v3 == v3) {
                triangles.remove(triangle)
                id = triangle.id
                break
            }
        }
        if(id != -1) {
            vertices[v1].adjacentTriangles.remove(id)
            vertices[v2].adjacentTriangles.remove(id)
            vertices[v3].adjacentTriangles.remove(id)

            if(makeVerticesUnneighbors(v1, v2)) removeEdge(v1, v2)
            if(makeVerticesUnneighbors(v1, v3)) removeEdge(v1, v3)
            if(makeVerticesUnneighbors(v2, v3)) removeEdge(v2, v3)
        }
    }

    fun splitTriangle(id: Int) {
        var triangle = Triangle(-1, -1, -1,-1)
        for(t in triangles) {
            if (t.id == id) {
                triangle = t
                break
            }
        }
        if(triangle.id == -1) return

        val centroid = addVertex(
            (vertices[triangle.v1].x + vertices[triangle.v2].x + vertices[triangle.v3].x) / 3,
            (vertices[triangle.v1].y + vertices[triangle.v2].y + vertices[triangle.v3].y) / 3,
            (vertices[triangle.v1].z + vertices[triangle.v2].z + vertices[triangle.v3].z) / 3
        )

        addTriangle(triangle.v1, triangle.v2, centroid)
        addTriangle(triangle.v2, triangle.v3, centroid)
        addTriangle(triangle.v3, triangle.v1, centroid)

        removeTriangle(triangle.v1, triangle.v2, triangle.v3)
    }

    fun split(line: String) : List<String> {
        return line.split(" ", "\t").filter{s -> s != ""}
    }

    fun loadOff(filename: String) {
        println("Initializing mesh using: $filename")
        val reader = File(filename).bufferedReader()

        // first line: OFF
        val header = reader.readLine()
        if(header != "OFF") {
            println("Unknown file format: $header")
            return
        }

        var line = reader.readLine()

        while(line.startsWith('#'))
            line = reader.readLine()

        // second line: #vertices #faces #edges
        val counters = split(line).map(String::toInt)

        for(i in 0 until counters[0]) {
            val coords = split(reader.readLine()).map(String::toFloat)
            addVertex(coords[0], coords[1], coords[2])
        }

        for(i in 0 until counters[1]) {
            val triangle = split(reader.readLine()).map(String::toInt)
            addTriangle(triangle[1], triangle[2], triangle[3])
        }

        println("Init completed. Mesh has ${triangles.size} triangles, ${vertices.size} vertices, ${edges.size} edges.")
    }

    fun exportOff(filename: String) {
        println("Exporting mesh to $filename")
        File(filename).printWriter().use { off ->
            off.println("OFF")
            off.println("${vertices.size} ${triangles.size} ${edges.size}")
            for (vertex in vertices) {
                off.println("${vertex.x.format()} ${vertex.y.format()} ${vertex.z.format()}")
            }
            for (triangle in triangles) {
                off.println("3 ${triangle.v1} ${triangle.v2} ${triangle.v3}")
            }
        }
    }

    private fun Float.format() = "%.10f".format(this)
}
