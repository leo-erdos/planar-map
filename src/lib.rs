use std::collections::LinkedList;

use rand::seq::SliceRandom;

struct Permutation {
    data: Vec<u64>,
}

impl Permutation {
    pub fn with_capacity(num_elements: usize) -> Self {
        Self {
            data: Vec::with_capacity(num_elements),
        }
    }

    pub fn cycle_push(&mut self, cycle_element: u64, new_element: u64) {
        let post = self.data[cycle_element as usize];
        self.data[cycle_element as usize] = new_element;
        self.data[new_element as usize] = post;
    }

    pub fn new_element(&mut self) -> u64 {
        let new_elt = self.data.len() as u64;
        self.data.push(new_elt);
        new_elt
    }

    pub fn transpose(&mut self, a: u64, b: u64) {
        let tmp = self.data[a as usize];
        self.data[a as usize] = self.data[b as usize];
        self.data[b as usize] = tmp;
    }
}

/// A planar map is an embedding of a graph into the sphere, up to orientation-preserving homeomorphisms.
pub struct PlanarMap {
    rotation_system: Permutation,
    edge_bindings: Permutation,
}

impl PlanarMap {
    /// Generates a random rooted planar tree with `n` nodes, using `rng` as its random number generator.
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
        let mut min = 0u64;
        let mut min_index = 0;
        let mut acc = 0u64;
        for (i, &step) in steps.iter().enumerate() {
            acc += step as u64;
            if acc < min {
                min = acc;
                min_index = i;
            }
        }
        // Build the tree
        let mut rotation_system = Permutation::with_capacity(2 * n);
        let mut edge_bindings = Permutation::with_capacity(2 * n);
        let mut vertex_stack: Vec<u64> = Vec::with_capacity(n);
        let mut acc = 0u64;
        for _i in 0..steps.len() {
            let i = _i + min_index % steps.len(); // Start at minimum and wrap around
            // Create new half edge
            let new_vertex = rotation_system.new_element();
            let new_half_edge = edge_bindings.new_element();
            debug_assert_eq!(new_vertex, new_half_edge);
            if acc < vertex_stack.len() as u64 {
                // Bind half edges
                let downward_half_edge = vertex_stack.pop().unwrap();
                let upward_half_edge = vertex_stack.pop().unwrap();
                edge_bindings.transpose(upward_half_edge, downward_half_edge);
                // Add to new half edge to vertex rotation system
                rotation_system.cycle_push(upward_half_edge, new_half_edge);
            }
            vertex_stack.push(new_vertex); // Add dangling half edge
            acc += steps[i] as u64;
        }
        Self {
            rotation_system: rotation_system,
            edge_bindings: edge_bindings,
        }
    }
}
