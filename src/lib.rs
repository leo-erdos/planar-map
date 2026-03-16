use bimap::BiMap;
use rand::seq::SliceRandom;
use std::{
    collections::{HashMap, HashSet, LinkedList, VecDeque},
    fmt::Debug,
    fs::File,
    io::BufWriter,
    io::Write,
};

#[derive(Debug, Clone)]
pub struct Permutation {
    data: Vec<u64>,
}

impl Permutation {
    pub fn data(&self) -> &Vec<u64> {
        &self.data
    }
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

    pub fn cycle_push(&mut self, cycle_element: u64, new_element: u64) -> &mut Self {
        let post = self.data[cycle_element as usize];
        self.data[cycle_element as usize] = new_element;
        self.data[new_element as usize] = post;
        self
    }

    pub fn cycles(&self) -> CyclesIterator {
        CyclesIterator::new(self)
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

    pub fn transpose(&mut self, a: u64, b: u64) -> &mut Self {
        let tmp = self.data[a as usize];
        self.data[a as usize] = self.data[b as usize];
        self.data[b as usize] = tmp;
        self
    }

    pub fn num_elements(&self) -> u64 {
        self.data.len() as u64
    }

    pub fn at(&self, element: u64) -> u64 {
        self.data[element as usize]
    }

    pub fn cycle_map(&self) -> Vec<u64> {
        let mut visited = vec![false; self.num_elements() as usize];
        let mut result = vec![0; self.num_elements() as usize];
        for i in 0..self.num_elements() {
            if visited[i as usize] {
                continue;
            }
            let mut v = i;
            loop {
                result[v as usize] = i;
                visited[v as usize] = true;
                v = self.at(v);
                if v == i {
                    break;
                }
            }
        }
        result
    }

    pub fn cycle_iter(&self, element: u64) -> CycleIterator {
        CycleIterator::new(self, element)
    }
}

pub struct CycleIterator<'a> {
    permutation: &'a Permutation,
    start: u64,
    current: u64,
    cycled: bool,
}

impl<'a> CycleIterator<'a> {
    pub fn new(permutation: &'a Permutation, start: u64) -> CycleIterator {
        CycleIterator {
            permutation: permutation,
            start: start,
            current: start,
            cycled: false,
        }
    }
}

impl<'a> Iterator for CycleIterator<'a> {
    type Item = u64;
    fn next(&mut self) -> Option<Self::Item> {
        match self.current == self.start && self.cycled {
            true => None,
            false => {
                self.cycled = true;
                let tmp = self.current;
                self.current = self.permutation.at(self.current);
                Some(tmp)
            }
        }
    }
}

pub struct CyclesIterator<'a> {
    permutation: &'a Permutation,
    visited: Vec<bool>,
    current: u64,
}

impl<'a> CyclesIterator<'a> {
    pub fn new(permutation: &'a Permutation) -> Self {
        Self {
            permutation: permutation,
            visited: vec![false; permutation.num_elements() as usize],
            current: 0,
        }
    }
}

impl<'a> Iterator for CyclesIterator<'a> {
    type Item = CycleIterator<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current == self.permutation.num_elements() {
            return None;
        }
        for cycle_element in self.permutation.cycle_iter(self.current) {
            self.visited[cycle_element as usize] = true;
        }
        let result = self.current;
        while self.current < self.permutation.num_elements() && self.visited[self.current as usize]
        {
            self.current += 1;
        }
        Some(self.permutation.cycle_iter(result))
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

    pub fn edges(&self) -> CyclesIterator {
        self.edge_bindings.cycles()
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

    pub fn vertices(&self) -> CyclesIterator {
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

// Helper for `tree_to_quadrangulation`. Used to connect all corners with label k+1
// to the following corner with label k
fn connect_corners_from_all(
    map: &mut PlanarMap,
    dst: u64,
    srcs: &Vec<u64>,
    inverse_rotation_system: &mut Permutation,
) {
    for &predecessor in srcs.iter().rev() {
        inverse_rotation_system.data[dst as usize] =
            connect_corners(map, dst, predecessor, inverse_rotation_system);
    }
}

/// Runs the CVS bijection on `tree`, using `labels`. Returns a quadrangulation.
pub fn tree_to_quadrangulation(tree: &mut PlanarMap, labels: &Vec<u64>) {
    assert_eq!(labels.len() as u64, tree.num_half_edges());
    let new_root = tree.add_half_edge();
    let mut inverse_rotation_system = tree.rotation_system().inverse(); // Used to calculate corners
    let mut label_stacks: Vec<Vec<u64>> =
        vec![Vec::new(); (tree.num_half_edges() / 2 + 3) as usize]; // Record of previously seen half-edges and their labels
    for i in 0..tree.num_half_edges() - 1 {
        let label = labels[i as usize] as usize;
        label_stacks[label].push(i);
        connect_corners_from_all(
            tree,
            i,
            &label_stacks[label + 1],
            &mut inverse_rotation_system,
        );
        label_stacks[label + 1].clear();
    }
    connect_corners_from_all(tree, 0, &label_stacks[2], &mut inverse_rotation_system);
    // Handle all vertices with label 1
    connect_corners_from_all(
        tree,
        new_root,
        &label_stacks[1],
        &mut inverse_rotation_system,
    );
    // Delete old edges from `tree`
    tree.rotation_system_mut().drop_prefix(new_root + 1);
    tree.edge_bindings_mut().drop_prefix(new_root + 1);
}

pub fn quadrangulation_to_triangulation(quad: &mut PlanarMap) -> &mut PlanarMap {
    let vertex_map = quad.rotation_system().cycle_map();
    let mut inverse_rotation_system = quad.rotation_system().inverse();
    let mut seen_vertices = HashSet::new();
    let mut face_vec = Vec::with_capacity(4);
    let face_walk = quad.rotation_system().compose(quad.edge_bindings());
    for face in face_walk.cycles() {
        face_vec.clear();
        let mut repeated = None;
        for (i, half_edge) in face.enumerate() {
            let vertex = vertex_map[half_edge as usize];
            if seen_vertices.contains(&vertex) {
                repeated = Some(i);
            }
            seen_vertices.insert(vertex);
            face_vec.push(half_edge);
        }
        seen_vertices.clear();
        match repeated {
            None => connect_corners(quad, face_vec[0], face_vec[2], &mut inverse_rotation_system),
            Some(i) => connect_corners(
                quad,
                face_vec[(i - 1) % 4],
                face_vec[(i + 1) % 4],
                &mut inverse_rotation_system,
            ),
        };
    }
    quad
}

pub fn is_simple(map: &mut PlanarMap) -> bool {
    use std::cmp::{max, min};
    let vertex_map = map.rotation_system().cycle_map();
    let mut edges: HashMap<(u64, u64), Vec<u64>> = HashMap::new();
    for half_edge in 0..map.num_half_edges() {
        let opposite_half_edge = map.edge_bindings().at(half_edge);
        if half_edge > opposite_half_edge {
            continue;
        }
        let vertex1 = vertex_map[half_edge as usize];
        let vertex2 = vertex_map[opposite_half_edge as usize];
        let abstract_edge = (min(vertex1, vertex2), max(vertex1, vertex2));
        match edges.get_mut(&abstract_edge) {
            None => {
                edges.insert(abstract_edge, vec![half_edge]);
            }
            Some(es) => {
                es.push(half_edge);
            }
        }
    }
    for edges_between in edges.values() {
        if edges_between.len() != 1 {
            return false;
        }
    }
    return true;
}

pub fn simplify_triangulation(triangulation: &mut PlanarMap) -> Vec<(u64, u64)> {
    use std::cmp::{max, min};
    let vertex_map = triangulation.rotation_system().cycle_map();
    let mut edges: HashMap<(u64, u64), Vec<u64>> = HashMap::new();
    let mut true_edges = Vec::new();
    for half_edge in 0..triangulation.num_half_edges() {
        let opposite_half_edge = triangulation.edge_bindings().at(half_edge);
        if half_edge > opposite_half_edge {
            continue;
        }
        let vertex1 = vertex_map[half_edge as usize];
        let vertex2 = vertex_map[opposite_half_edge as usize];
        let abstract_edge = (min(vertex1, vertex2), max(vertex1, vertex2));
        match edges.get_mut(&abstract_edge) {
            None => {
                edges.insert(abstract_edge, vec![half_edge]);
            }
            Some(es) => {
                es.push(half_edge);
            }
        }
    }
    for edges_between in edges.values() {
        if edges_between.len() == 1 {
            continue;
        }
        for &half_edge in edges_between {
            // Add a new vertex in the middle of the edge
            let new1 = triangulation.add_half_edge();
            let new2 = triangulation.add_half_edge();
            let opposite_half_edge = triangulation.edge_bindings().at(half_edge);
            triangulation.rotation_system_mut().transpose(new1, new2);
            triangulation
                .edge_bindings_mut()
                .transpose(half_edge, opposite_half_edge)
                .transpose(half_edge, new1)
                .transpose(opposite_half_edge, new2);
            true_edges.push((half_edge, opposite_half_edge));
            // Connect the new vertex to the non-adjacent vertex in both adjacent faces
            let new3 = triangulation.add_half_edge();
            let new4 = triangulation.add_half_edge();
            let new5 = triangulation.add_half_edge();
            let new6 = triangulation.add_half_edge();
            let c1 = triangulation
                .edge_bindings()
                .at(triangulation.rotation_system().at(half_edge));
            let c2 = triangulation
                .edge_bindings()
                .at(triangulation.rotation_system().at(opposite_half_edge));
            triangulation
                .rotation_system_mut()
                .cycle_push(new1, new4)
                .cycle_push(new2, new3)
                .cycle_push(c1, new5)
                .cycle_push(c2, new6);
            triangulation
                .edge_bindings_mut()
                .transpose(new3, new5)
                .transpose(new4, new6);
        }
    }
    true_edges
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum HalfEdgeColor {
    None,
    Red,
    Blue,
    Green,
}

#[derive(Clone, Copy, Debug)]
struct Corner {
    start: u64,
    end: u64,
}

enum OuterCycleElement {
    Corner(Corner, u8),
    Redirect(usize, usize),
}

struct OuterCycle {
    elements: VecDeque<OuterCycleElement>,
    index: usize,
}

impl OuterCycle {
    pub fn new(initial_corners: &[Corner]) -> Self {
        let mut elements = VecDeque::with_capacity(initial_corners.len());
        for &corner in initial_corners {
            elements.push_back(OuterCycleElement::Corner(corner, 0));
        }
        Self {
            elements: elements,
            index: 2,
        }
    }

    fn get(&mut self) -> &mut Corner {
        match self.elements.get_mut(self.index) {
            Some(OuterCycleElement::Corner(corner, _)) => corner,
            _ => panic!("unreachable"),
        }
    }

    pub fn next(&mut self) -> &mut Corner {
        self.index += 1;
        if self.index >= self.elements.len() {
            self.index = 2;
        }
        while let OuterCycleElement::Redirect(new_index, _) = self.elements[self.index] {
            self.index = new_index;
        }
        self.get()
    }

    pub fn previous(&mut self) -> &mut Corner {
        self.index -= 1;
        if self.index < 2 {
            self.index = self.elements.len() - 1;
        }
        while let OuterCycleElement::Redirect(_, new_index) = self.elements[self.index] {
            self.index = new_index;
        }
        self.get()
    }

    /*pub fn next_nonchordal(&mut self) -> &mut Corner {
        match self.elements[self.index] {
            OuterCycleElement::Redirect(_) => {
                self.next();
            }
            _ => (),
        };
        while match self.elements[self.index] {
            OuterCycleElement::Corner(_, chords) => chords != 0,
            _ => true,
        } {
            self.next();
        }
        self.get()
    }*/

    /*fn expand_cycle(&mut self, new_corners: &[Corner], chordalities: &[u8]) -> &mut Self {
        self.
        for i in new_corners.len() {

        }
    }*/
}

fn has_chords(
    triangulation: &PlanarMap,
    corner: (u64, u64),
    outer_cycle: &BiMap<u64, u64>,
    vertex_map: &Vec<u64>,
    cycle_vertices: &HashSet<u64>,
) -> bool {
    if outer_cycle.len() <= 3 {
        return false;
    }
    let mut h = triangulation.next_rotation(corner.0);
    while h != corner.1 {
        let neighbor = vertex_map[triangulation.edge_bindings().at(h) as usize];
        if cycle_vertices.contains(&neighbor) {
            return true;
        }
        h = triangulation.next_rotation(h);
    }
    false
}

pub fn schnyder_woods(triangulation: &PlanarMap) -> Vec<HalfEdgeColor> {
    use HalfEdgeColor::*;
    let mut colors = vec![None; triangulation.num_half_edges() as usize];
    let vertex_map = triangulation.rotation_system().cycle_map();
    let inverse_rotation_system = triangulation.rotation_system().inverse();
    let face_walk = triangulation
        .rotation_system()
        .compose(triangulation.edge_bindings());
    let mut outer_cycle: BiMap<u64, u64> = BiMap::new();
    let g = face_walk.at(0);
    let b = face_walk.at(g);
    let red_vertex = vertex_map[0];
    colors[triangulation.edge_bindings().at(0) as usize] = HalfEdgeColor::Red;
    colors[g as usize] = HalfEdgeColor::Blue;
    colors[triangulation.edge_bindings().at(g) as usize] = HalfEdgeColor::Green;
    colors[b as usize] = HalfEdgeColor::Red;
    let green_vertex = vertex_map[g as usize];
    let blue_vertex = vertex_map[b as usize];
    let mut cycle_vertices: HashSet<u64> = HashSet::from([red_vertex, green_vertex, blue_vertex]);
    for half_edge in face_walk.cycle_iter(0) {
        outer_cycle.insert(half_edge, inverse_rotation_system.at(half_edge));
    }
    let mut corner = 0;
    while outer_cycle.len() > 2 {
        let &corner_end = outer_cycle.get_by_left(&corner).unwrap();
        if vertex_map[corner as usize] == blue_vertex || vertex_map[corner as usize] == green_vertex
        {
            corner = triangulation.edge_bindings().at(corner_end);
            continue;
        }
        if has_chords(
            triangulation,
            (corner, corner_end),
            &outer_cycle,
            &vertex_map,
            &cycle_vertices,
        ) {
            corner = triangulation.edge_bindings().at(corner_end);
            continue;
        }
        colors[corner as usize] = Green;
        colors[corner_end as usize] = Blue;
        let mut to_neighbor = triangulation.rotation_system().at(corner);
        let mut previous_neighbor = triangulation.edge_bindings().at(corner);
        outer_cycle.remove_by_left(&corner);
        cycle_vertices.remove(&vertex_map[corner as usize]);
        while to_neighbor != corner_end {
            let neighbor = triangulation.edge_bindings().at(to_neighbor);
            colors[neighbor as usize] = Red;
            outer_cycle.insert(triangulation.rotation_system().at(neighbor), neighbor);
            cycle_vertices.insert(vertex_map[neighbor as usize]);
            let (previous_neighbor_start, _) =
                outer_cycle.remove_by_right(&previous_neighbor).unwrap();
            outer_cycle.insert(
                previous_neighbor_start,
                inverse_rotation_system.at(previous_neighbor),
            );
            previous_neighbor = neighbor;
            to_neighbor = triangulation.rotation_system().at(to_neighbor);
        }
        let last_neighbor = triangulation.edge_bindings().at(corner_end);
        let (_, last_neighbor_end) = outer_cycle.remove_by_left(&last_neighbor).unwrap();
        outer_cycle.insert(
            triangulation.rotation_system().at(last_neighbor),
            last_neighbor_end,
        );
        let (previous_neighbor_start, _) = outer_cycle.remove_by_right(&previous_neighbor).unwrap();
        outer_cycle.insert(
            previous_neighbor_start,
            inverse_rotation_system.at(previous_neighbor),
        );
        corner = triangulation.rotation_system().at(last_neighbor);
    }
    colors
}

fn colored_parents<'a>(
    triangulation: &'a PlanarMap,
    colors: &'a Vec<HalfEdgeColor>,
    vertex: u64,
    color: HalfEdgeColor,
    vertex_map: &'a Vec<u64>,
) -> impl Iterator<Item = u64> + 'a {
    triangulation
        .rotation_system()
        .cycle_iter(vertex)
        .map(|half_edge| (triangulation.edge_bindings().at((half_edge))))
        .filter(move |&opp_half_edge| (colors[opp_half_edge as usize]) == (color))
        .map(|opp_half_edge| vertex_map[opp_half_edge as usize])
}

fn count_ancestors(
    triangulation: &PlanarMap,
    colors: &Vec<HalfEdgeColor>,
    vertex: u64,
    color: HalfEdgeColor,
    counts: &mut HashMap<u64, u64>,
    vertex_map: &Vec<u64>,
) -> u64 {
    let vertex = vertex_map[vertex as usize];
    if counts.contains_key(&vertex) {
        return counts[&vertex];
    }
    let mut count = 0;
    for parent in colored_parents(triangulation, colors, vertex, color, vertex_map) {
        count += count_ancestors(triangulation, colors, parent, color, counts, vertex_map) + 1;
    }
    counts.insert(vertex, count);
    count
}

fn colored_child(
    triangulation: &PlanarMap,
    colors: &Vec<HalfEdgeColor>,
    vertex: u64,
    color: HalfEdgeColor,
    vertex_map: &Vec<u64>,
) -> Option<u64> {
    for half_edge in triangulation.rotation_system().cycle_iter(vertex) {
        if colors[half_edge as usize] == color {
            return Some(vertex_map[triangulation.edge_bindings.at(half_edge) as usize]);
        }
    }
    None
}

fn count_descendants(
    triangulation: &PlanarMap,
    colors: &Vec<HalfEdgeColor>,
    vertex: u64,
    color: HalfEdgeColor,
    counts: &mut HashMap<u64, u64>,
    vertex_map: &Vec<u64>,
) -> u64 {
    let vertex = vertex_map[vertex as usize];
    if counts.contains_key(&vertex) {
        return counts[&vertex];
    }
    let count = match colored_child(triangulation, colors, vertex, color, vertex_map) {
        None => 0,
        Some(child) => {
            count_descendants(triangulation, colors, child, color, counts, vertex_map) + 1
        }
    };
    counts.insert(vertex, count);
    count
}

fn count_ancestors_on_path(
    triangulation: &PlanarMap,
    colors: &Vec<HalfEdgeColor>,
    vertex: u64,
    ancestor_color: HalfEdgeColor,
    path_color: HalfEdgeColor,
    counts: &mut HashMap<u64, u64>,
    ancestor_counts: &mut HashMap<u64, u64>,
    vertex_map: &Vec<u64>,
) -> u64 {
    let vertex = vertex_map[vertex as usize];
    if counts.contains_key(&vertex) {
        return counts[&vertex];
    }
    let count = count_ancestors(
        triangulation,
        colors,
        vertex,
        ancestor_color,
        ancestor_counts,
        vertex_map,
    ) + match colored_child(triangulation, colors, vertex, path_color, vertex_map) {
        None => 0,
        Some(child) => count_ancestors_on_path(
            triangulation,
            colors,
            child,
            ancestor_color,
            path_color,
            counts,
            ancestor_counts,
            vertex_map,
        ),
    };
    counts.insert(vertex, count);
    count
}

pub fn schnyder_embedding(
    triangulation: &PlanarMap,
    colors: &Vec<HalfEdgeColor>,
) -> HashMap<u64, (i64, i64)> {
    let vertex_map = triangulation.rotation_system().cycle_map();
    let mut green_ancestor_counts = HashMap::new();
    let mut blue_ancestor_counts = HashMap::new();
    let green_sink = triangulation
        .rotation_system()
        .at(triangulation.edge_bindings().at(0));
    let blue_sink = triangulation
        .rotation_system()
        .at(triangulation.edge_bindings().at(green_sink));
    count_ancestors(
        triangulation,
        colors,
        green_sink,
        HalfEdgeColor::Green,
        &mut green_ancestor_counts,
        &vertex_map,
    );
    count_ancestors(
        triangulation,
        colors,
        blue_sink,
        HalfEdgeColor::Blue,
        &mut blue_ancestor_counts,
        &vertex_map,
    );
    let mut blue_descendant_count = HashMap::new();
    let mut green_descendant_count = HashMap::new();
    let mut red_descendant_count = HashMap::new();
    let mut red_path_green_ancestor_counts = HashMap::new();
    let mut blue_path_green_ancestor_counts = HashMap::new();
    let mut red_path_blue_ancestor_counts = HashMap::new();
    let mut green_path_blue_ancestor_counts = HashMap::new();
    for x in 0..triangulation.num_half_edges() {
        count_descendants(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Red,
            &mut red_descendant_count,
            &vertex_map,
        );
        count_descendants(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Green,
            &mut green_descendant_count,
            &vertex_map,
        );
        count_descendants(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Blue,
            &mut blue_descendant_count,
            &vertex_map,
        );
        count_ancestors_on_path(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Green,
            HalfEdgeColor::Red,
            &mut red_path_green_ancestor_counts,
            &mut green_ancestor_counts,
            &vertex_map,
        );
        count_ancestors_on_path(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Green,
            HalfEdgeColor::Blue,
            &mut blue_path_green_ancestor_counts,
            &mut green_ancestor_counts,
            &vertex_map,
        );
        count_ancestors_on_path(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Blue,
            HalfEdgeColor::Red,
            &mut red_path_blue_ancestor_counts,
            &mut blue_ancestor_counts,
            &vertex_map,
        );
        count_ancestors_on_path(
            triangulation,
            colors,
            x,
            HalfEdgeColor::Blue,
            HalfEdgeColor::Green,
            &mut green_path_blue_ancestor_counts,
            &mut blue_ancestor_counts,
            &vertex_map,
        );
    }
    let vertices = green_ancestor_counts[&vertex_map[green_sink as usize]] as i64;
    let mut positions: HashMap<u64, (i64, i64)> = HashMap::with_capacity(vertices as usize);
    positions.insert(0, (0, vertices));
    positions.insert(vertex_map[green_sink as usize], (vertices, -vertices));
    positions.insert(vertex_map[blue_sink as usize], (-vertices, -vertices));
    for half_edge in 0..triangulation.num_half_edges() {
        let vertex = vertex_map[half_edge as usize];
        if positions.contains_key(&vertex) {
            continue;
        }
        let green_coordinate = (blue_path_green_ancestor_counts[&vertex]
            + red_path_green_ancestor_counts[&vertex]
            + red_descendant_count[&vertex]) as i64
            - (green_ancestor_counts[&vertex] as i64);
        let blue_coordinate = (red_path_blue_ancestor_counts[&vertex]
            + green_path_blue_ancestor_counts[&vertex]
            + green_descendant_count[&vertex]) as i64
            - (blue_ancestor_counts[&vertex] as i64);
        let red_coordinate = vertices - green_coordinate - blue_coordinate;
        positions.insert(
            vertex,
            (
                green_coordinate - blue_coordinate,
                red_coordinate - green_coordinate - blue_coordinate,
            ),
        );
    }
    positions
}

pub fn initial_3d_embedding<R: rand::Rng>(
    positions: &HashMap<u64, (i64, i64)>,
    vertices: u64,
    rng: &mut R,
    fixed: &HashSet<u64>,
) -> HashMap<u64, (f64, f64, f64)> {
    let vertices = vertices as f64;
    let mut result = HashMap::new();
    for (&vertex, &position) in positions.iter() {
        if fixed.contains(&vertex) {
            result.insert(
                vertex,
                (
                    (position.0 as f64) / vertices,
                    (position.1 as f64) / vertices,
                    0.0,
                ),
            );
            continue;
        }
        let x = (position.0 as f64) / vertices;
        let y = (position.1 as f64) / vertices;
        let z = (1.0 - (x * x + y * y).sqrt()) / 2.0 + rng.random::<f64>() / 4.0;
        result.insert(vertex, (x, y, z));
    }
    result
}

pub fn save_obj(triangulation: &PlanarMap, positions: &HashMap<u64, (f64, f64, f64)>, path: &str) {
    let file = File::create(path).unwrap();
    let mut w = BufWriter::new(file);
    let vertex_index = vertex_index_map(triangulation);
    for i in 0..vertex_index.len() {
        let vertex = vertex_index.get_by_right(&(i as u64)).unwrap();
        let (x, y, z) = positions[&vertex];
        writeln!(w, "v {x} {y} {z}").unwrap();
    }
    let vertex_map = triangulation.rotation_system().cycle_map();
    let face_walk = triangulation
        .rotation_system()
        .compose(triangulation.edge_bindings());
    for face in face_walk.cycles() {
        write!(w, "f").unwrap();
        for half_edge in face {
            let vertex = vertex_map[half_edge as usize];
            let &line = vertex_index.get_by_left(&vertex).unwrap();
            let line = line + 1;
            write!(w, " {line}").unwrap();
        }
        write!(w, "\n").unwrap();
    }
}

pub fn distance_matrix(map: &PlanarMap) -> Vec<u64> {
    let n_vertices = map.rotation_system().cycles().count();
    let vertex_map = map.rotation_system().cycle_map();
    let mut vertex_index_map = HashMap::with_capacity(n_vertices);
    let mut result = vec![0; n_vertices * n_vertices];
    let mut visiting = Vec::with_capacity(n_vertices);
    let mut i = 0;
    for mut vertex in map.rotation_system().cycles() {
        let representant = vertex.next().unwrap();
        visiting.push((representant, representant));
        vertex_index_map.insert(representant, i);
        i += 1;
    }
    let mut new_visiting;
    let n_vertices = n_vertices as u64;
    while !visiting.is_empty() {
        new_visiting = Vec::with_capacity(n_vertices as usize);
        for &(source, destination) in visiting.iter() {
            let dist = result[(vertex_index_map[&source] * n_vertices
                + vertex_index_map[&destination]) as usize];
            for half_edge in map.rotation_system().cycle_iter(destination) {
                let neighbor = vertex_map[map.edge_bindings().at(half_edge) as usize];
                if result[(vertex_index_map[&source] * n_vertices + vertex_index_map[&neighbor])
                    as usize]
                    > dist + 1
                {
                    result[(vertex_index_map[&source] * n_vertices + vertex_index_map[&neighbor])
                        as usize] = dist + 1;
                }
                if result[(vertex_index_map[&source] * n_vertices + vertex_index_map[&neighbor])
                    as usize]
                    == 0
                    && neighbor != source
                {
                    result[(vertex_index_map[&source] * n_vertices + vertex_index_map[&neighbor])
                        as usize] = dist + 1;
                    new_visiting.push((source, neighbor));
                }
            }
        }
        visiting = new_visiting;
    }
    result
}

pub fn vertex_index_map(map: &PlanarMap) -> BiMap<u64, u64> {
    let mut result = BiMap::new();
    let mut i = 0;
    for mut vertex_iter in map.rotation_system().cycles() {
        result.insert(vertex_iter.next().unwrap(), i);
        i += 1;
    }
    result
}

pub fn augmented_distance_matrix(
    map: &PlanarMap,
    distance_matrix: &Vec<u64>,
    n_vertices_no_augment: u64,
    n_half_edges_no_augment: u64,
) -> Vec<u64> {
    let n_vertices = map.rotation_system().cycles().count();
    let mut result = vec![0; n_vertices * n_vertices];
    for i in 0..distance_matrix.len() {
        let row = i / n_vertices_no_augment as usize;
        let col = i % n_vertices_no_augment as usize;
        result[row * n_vertices + col] = 2 * distance_matrix[i];
    }
    let vertex_map = map.rotation_system().cycle_map();
    let vertex_index = vertex_index_map(map);
    let mut index_neighbors_map: Vec<Vec<u64>> =
        Vec::with_capacity(n_vertices - n_vertices_no_augment as usize);
    let mut i = 0;
    for mut vertex_iter in map.rotation_system().cycles() {
        let representant = vertex_iter.next().unwrap();
        if *vertex_index.get_by_left(&representant).unwrap() >= n_vertices_no_augment {
            let our_neighbors: Vec<u64> = map
                .rotation_system()
                .cycle_iter(representant)
                .map(|x| map.edge_bindings().at(x))
                .filter(|x| *x < n_half_edges_no_augment)
                .map(|x| vertex_map[x as usize])
                .collect();
            index_neighbors_map.push(our_neighbors.clone());
            for k in 0..i {
                let their_neighbors = &index_neighbors_map[k];
                let mut min_dist =
                    distance_matrix[(*vertex_index.get_by_left(&our_neighbors[0]).unwrap()
                        * n_vertices_no_augment as u64
                        + *vertex_index.get_by_left(&their_neighbors[0]).unwrap())
                        as usize];
                for n1 in our_neighbors.iter() {
                    for n2 in their_neighbors.iter() {
                        let dist = distance_matrix[((vertex_index.get_by_left(n1)).unwrap()
                            * n_vertices_no_augment as u64
                            + vertex_index.get_by_left(n2).unwrap())
                            as usize];
                        if dist < min_dist {
                            min_dist = dist;
                        }
                    }
                }
                result[((n_vertices_no_augment + i as u64) * n_vertices as u64
                    + (n_vertices_no_augment + k as u64)) as usize] = 2 * (min_dist + 1);
                result[((n_vertices_no_augment + k as u64) * n_vertices as u64
                    + (n_vertices_no_augment + i as u64)) as usize] = 2 * (min_dist + 1);
            }
            for x in 0..n_vertices_no_augment {
                let min_dist =
                    (distance_matrix[(vertex_index.get_by_left(&our_neighbors[0]).unwrap()
                        * n_vertices_no_augment as u64
                        + x) as usize])
                        .min(
                            distance_matrix[(vertex_index.get_by_left(&our_neighbors[1]).unwrap()
                                * n_vertices_no_augment as u64
                                + x) as usize],
                        );
                result[((n_vertices_no_augment + i as u64) * n_vertices as u64 + x) as usize] =
                    2 * min_dist + 1;
                result[(x * n_vertices as u64 + n_vertices_no_augment + i as u64) as usize] =
                    2 * min_dist + 1;
            }
            i += 1;
        }
    }
    result
}

pub fn save_quad_obj(map: &PlanarMap, diag_start: u64, true_edges: &Vec<(u64, u64)>, path: &str) {
    let vertex_index = vertex_index_map(map);
    let vertex_map = map.rotation_system().cycle_map();
    let mut file = File::create(path).unwrap();
    for edge_iter in map.edges() {
        let edge: Vec<u64> = edge_iter.collect();
        let h0 = edge[0];
        let h1 = edge[1];
        if h0 < diag_start && h1 < diag_start {
            write!(
                file,
                "{} {}\n",
                vertex_index.get_by_left(&vertex_map[h0 as usize]).unwrap() + 1,
                vertex_index.get_by_left(&vertex_map[h1 as usize]).unwrap() + 1
            )
            .unwrap();
        }
    }
    for true_edge in true_edges {
        let h0 = true_edge.0;
        let h1 = true_edge.1;
        let h2 = map.edge_bindings().at(h0);
        let h3 = map.edge_bindings().at(h1);
        write!(
            file,
            "{} {} {} {}\n",
            vertex_index.get_by_left(&vertex_map[h0 as usize]).unwrap() + 1,
            vertex_index.get_by_left(&vertex_map[h2 as usize]).unwrap() + 1,
            vertex_index.get_by_left(&vertex_map[h3 as usize]).unwrap() + 1,
            vertex_index.get_by_left(&vertex_map[h1 as usize]).unwrap() + 1,
        )
        .unwrap();
    }
}
