use std::collections::LinkedList;

use rand::seq::SliceRandom;

#[derive(Debug, Clone)]
struct Permutation {
    data: Vec<u64>,
}

impl Permutation {
    pub fn with_capacity(num_elements: usize) -> Self {
        Self {
            data: Vec::with_capacity(num_elements),
        }
    }

    pub fn compose(&self, other: &Self) -> Self {
        assert_eq!(self.num_elements(), other.num_elements());
        let mut result = Self::with_capacity(self.num_elements() as usize);
        for i in 0..self.num_elements() {
            result.data.push(self.at(other.at(i)));
        }
        result
    }

    pub fn cycle_push(&mut self, cycle_element: u64, new_element: u64) {
        let post = self.data[cycle_element as usize];
        self.data[cycle_element as usize] = new_element;
        self.data[new_element as usize] = post;
    }

    pub fn cycles(&self) -> Vec<Vec<u64>> {
        let mut visited = vec![false; self.num_elements() as usize];
        let mut cycles = Vec::new();
        for i in 0..self.num_elements() {
            if visited[i as usize] {
                continue;
            }
            cycles.push(Vec::with_capacity(1));
            let current_cycle = cycles.last_mut().unwrap();
            let mut cycle_element = i;
            while !visited[cycle_element as usize] {
                current_cycle.push(cycle_element);
                visited[cycle_element as usize] = true;
                cycle_element = self.at(cycle_element);
            }
        }
        cycles
    }

    pub fn drop_prefix(&mut self, prefix_len: u64) {
        for i in prefix_len..self.num_elements() {
            while self.at(i) < prefix_len {
                self.data[i as usize] = self.at(self.at(i));
            }
        }
        self.data.drain(0..prefix_len as usize);
        for e in self.data.iter_mut() {
            *e -= prefix_len;
        }
    }

    pub fn inverse(&self) -> Self {
        let mut result_data = vec![0; self.num_elements() as usize];
        for i in 0..self.num_elements() {
            result_data[self.at(i) as usize] = i;
        }
        Self { data: result_data }
    }

    pub fn add_element(&mut self) -> u64 {
        let new_elt = self.data.len() as u64;
        self.data.push(new_elt);
        new_elt
    }

    pub fn transpose(&mut self, a: u64, b: u64) {
        let tmp = self.data[a as usize];
        self.data[a as usize] = self.data[b as usize];
        self.data[b as usize] = tmp;
    }

    pub fn num_elements(&self) -> u64 {
        self.data.len() as u64
    }

    pub fn at(&self, element: u64) -> u64 {
        self.data[element as usize]
    }
}

/// A planar map is an embedding of a graph into the sphere, up to orientation-preserving homeomorphisms.
#[derive(Clone, Debug)]
pub struct PlanarMap {
    rotation_system: Permutation,
    edge_bindings: Permutation,
}

impl PlanarMap {
    pub fn add_half_edge(&mut self) -> u64 {
        self.rotation_system.add_element();
        self.edge_bindings.add_element()
    }

    /// Returns the map's edge bindings as an immutable reference. Let e be an edge, composed of 2 half edges
    /// a and b. The edge bindings are a permutation that maps a -> b, and b -> a.
    pub fn edge_bindings(&self) -> &Permutation {
        &self.edge_bindings
    }

    /// Returns the map's edge bindings as an mutable reference. Let e be an edge, composed of 2 half edges
    /// a and b. The edge bindings are a permutation that maps a -> b, and b -> a.
    pub fn edge_bindings_mut(&mut self) -> &mut Permutation {
        &mut self.edge_bindings
    }

    pub fn edges(&self) -> Vec<Vec<u64>> {
        self.edge_bindings.cycles()
    }

    pub fn faces(&self) -> Vec<Vec<u64>> {
        self.rotation_system.compose(&self.edge_bindings).cycles()
    }

    /// Returns the half edge that follows `half_edge` in the rotation order, around the vertex
    /// incident to `half_edge`.
    pub fn next_rotation(&self, half_edge: u64) -> u64 {
        self.rotation_system.at(half_edge)
    }

    pub fn num_half_edges(&self) -> u64 {
        self.rotation_system.num_elements()
    }

    /// Generates a random rooted planar tree with `n` vertices, using `rng` as its random number generator.
    pub fn random_tree<R: rand::Rng>(n: usize, rng: &mut R) -> Self {
        // Generate a random sequence of upward (1) or downward (-1) steps with sum 0.
        let mut steps: Vec<i8> = Vec::with_capacity(2 * n);
        for _ in 0..n {
            steps.push(1);
        }
        for _ in 0..n {
            steps.push(-1);
        }
        steps.shuffle(rng);
        // Find the first minimum
        let mut min = 0i64;
        let mut min_index = 0;
        let mut acc = 0i64;
        for (i, &step) in steps.iter().enumerate() {
            if acc < min {
                min = acc;
                min_index = i;
            }
            acc += step as i64;
        }
        // Build the tree
        let mut rotation_system = Permutation::with_capacity(2 * n);
        let mut edge_bindings = Permutation::with_capacity(2 * n);
        let mut vertex_stack: Vec<u64> = Vec::with_capacity(n);
        let mut acc = 0i64;
        for _i in 0..steps.len() {
            let i = (_i + min_index) % steps.len(); // Start at minimum and wrap around
            // Create new half edge
            let new_vertex = rotation_system.add_element();
            let new_half_edge = edge_bindings.add_element();
            debug_assert_eq!(new_vertex, new_half_edge);
            if acc < vertex_stack.len() as i64 {
                // Bind half edges
                let downward_half_edge = vertex_stack.pop().unwrap();
                let upward_half_edge = vertex_stack.pop().unwrap();
                edge_bindings.transpose(upward_half_edge, downward_half_edge);
                // Add to new half edge to vertex rotation system
                rotation_system.cycle_push(upward_half_edge, new_half_edge);
            }
            vertex_stack.push(new_vertex); // Add dangling half edge
            acc += steps[i] as i64;
        }
        // Bind last dangling half-edges
        let downward_half_edge = vertex_stack.pop().unwrap();
        let upward_half_edge = vertex_stack.pop().unwrap();
        edge_bindings.transpose(upward_half_edge, downward_half_edge);
        Self {
            rotation_system: rotation_system,
            edge_bindings: edge_bindings,
        }
    }

    /// Returns the map's rotation system as an immutable reference. The rotation system is a permutation that maps
    /// any half-edge to the half-edge following in the rotation order around its incident vertex.
    pub fn rotation_system(&self) -> &Permutation {
        &self.rotation_system
    }

    /// Returns the map's rotation system as a mutable reference. The rotation system is a permutation that maps
    /// any half-edge to the half-edge following in the rotation order around its incident vertex.
    pub fn rotation_system_mut(&mut self) -> &mut Permutation {
        &mut self.rotation_system
    }

    pub fn vertices(&self) -> Vec<Vec<u64>> {
        self.rotation_system.cycles()
    }
}

// Helper for `labels`, sets the label for every half edge incident to `vertex`
fn populate_vertex_labels(tree: &PlanarMap, vertex: u64, new_label: u64, labels: &mut Vec<u64>) {
    let mut corner = vertex;
    loop {
        labels[corner as usize] = new_label;
        corner = tree.next_rotation(corner);
        if corner == vertex {
            break;
        }
    }
}

/// Computes vertex labels for the CVS bijection, using `rng` as its random number generator.
/// Returns a vector mapping half-edge numbers to their labels.
pub fn labels<R: rand::Rng>(tree: &PlanarMap, rng: &mut R) -> Vec<u64> {
    let mut result = vec![0; tree.num_half_edges() as usize]; // 0 means label not computed yet
    populate_vertex_labels(tree, 0, 1, &mut result);
    for i in 1..tree.num_half_edges() {
        if result[i as usize] == 0 {
            let label = match result[(i - 1) as usize] {
                1 => rng.random_range(1..=2),
                k => rng.random_range(k - 1..=k + 1),
            };
            populate_vertex_labels(tree, i, label, &mut result);
        }
    }
    result
}

// Helper for `tree_to_quadrangulation`. Connects 2 corners by a new edge
fn connect_corners(
    map: &mut PlanarMap,
    a: u64,
    b: u64,
    inverse_rotation_system: &Permutation,
) -> u64 {
    // Create half-edges
    let new_half_edge_a = map.add_half_edge();
    let new_half_edge_b = map.add_half_edge();
    // Add the half-edges in their respective corners
    map.rotation_system_mut()
        .cycle_push(inverse_rotation_system.at(a), new_half_edge_a);
    map.rotation_system_mut()
        .cycle_push(inverse_rotation_system.at(b), new_half_edge_b);
    // Bind half-edges
    map.edge_bindings_mut()
        .transpose(new_half_edge_a, new_half_edge_b);
    new_half_edge_a
}

/// Runs the CVS bijection on `tree`, using `labels`. Returns a quadrangulation.
pub fn tree_to_quadrangulation(tree: &PlanarMap, labels: &Vec<u64>) -> PlanarMap {
    assert_eq!(labels.len() as u64, tree.num_half_edges());
    let mut result = tree.clone();
    let new_root = result.add_half_edge();
    let mut inverse_rotation_system = result.rotation_system().inverse(); // Used to calculate corners
    let mut label_stacks: Vec<Vec<u64>> =
        vec![Vec::new(); (tree.num_half_edges() / 2 + 3) as usize]; // Record of previously seen half-edges and their labels
    for i in 0..tree.num_half_edges() {
        let label = labels[i as usize] as usize;
        label_stacks[label].push(i);
        let mut first_new_half_edge = None;
        for &predecessor in label_stacks[label + 1].iter() {
            // Add an edge towards the predecessor
            let new_half_edge =
                connect_corners(&mut result, i, predecessor, &inverse_rotation_system);
            if predecessor == *label_stacks[label + 1].first().unwrap() {
                first_new_half_edge = Some(new_half_edge);
            }
        }
        if let Some(x) = first_new_half_edge {
            inverse_rotation_system.data[i as usize] = x;
        }
        label_stacks[label + 1].clear();
    }
    let mut first_new_half_edge = None;
    // The initial half-edge may also be the successor of a half-edge
    for &predecessor in label_stacks[2].iter() {
        let new_half_edge = connect_corners(&mut result, 0, predecessor, &inverse_rotation_system);
        if predecessor == *label_stacks[2].first().unwrap() {
            first_new_half_edge = Some(new_half_edge);
        }
    }
    if let Some(x) = first_new_half_edge {
        inverse_rotation_system.data[0] = x;
    }
    //inverse_rotation_system = result.rotation_system().inverse();
    // Handle all vertices with label 1
    for &i in label_stacks[1].iter() {
        connect_corners(&mut result, new_root, i, &inverse_rotation_system);
    } // Delete old edges from `tree`
    result.rotation_system_mut().drop_prefix(new_root + 1);
    result.edge_bindings_mut().drop_prefix(new_root + 1);
    result
}
