aabb-tree
=========

A d-dimensional aabb-tree implementation for MATLAB.

AABBTREE offers d-dimensional aabb-tree construction and search for collections of spatial objects. These trees are useful when seeking to implement efficient spatial queries -- determining intersections between collections of objects.

Given a collection of objects, an aabb-tree partitions the bounding-boxes associated with the elements in the collection into a (binary) "tree" -- a hierarchy of rectangular "nodes" that each store a subset of the collection. In contrast to other geometric tree types (quadtrees, kd-trees, etc), aabb-trees are applicable to collections of general objects, rather than just points.

AABBTREE is expected to support additional search queries in the future.

See AABBDEMO, MAKETREE for additional details.
