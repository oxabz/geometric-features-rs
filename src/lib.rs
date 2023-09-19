/*!
Rust Library Dedicated to compute features of an object based on it's polygon coordinates
*/
pub mod convex_hull;
#[cfg(feature = "ellipses")]
pub mod related_ellipse;
mod utils;
use std::f64::consts;

pub type Point = (f64, f64);
pub type Polygon<'a> = &'a [Point];
pub type PolygonOwned = Vec<Point>;

/**
Compute the area of a polygon

```rust
use geometric_features::area;
let polygon: Vec<(f64,f64)> = vec![
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.0),
    (0.0, 1.0),
];
assert_eq!(area(&polygon), 1.0);
```
 */
pub fn area(polygon: Polygon) -> f64 {
    let mut area = 0.0;
    let mut j = polygon.len() - 1;
    for i in 0..polygon.len() {
        area = area + (polygon[j].0 + polygon[i].0) * (polygon[j].1 - polygon[i].1);
        j = i;
    }
    area.abs() / 2.0
}

/**
Compute the perimeter of a polygon

```rust
use geometric_features::perimeter;
let polygon: Vec<(f64,f64)> = vec![
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.0),
    (0.0, 1.0),
];
assert_eq!(perimeter(&polygon), 4.0);
```
 */
pub fn perimeter(polygon: Polygon) -> f64 {
    let mut perimeter = 0.0;
    let mut j = polygon.len() - 1;
    for i in 0..polygon.len() {
        perimeter = perimeter
            + ((polygon[j].0 - polygon[i].0).powi(2) + (polygon[j].1 - polygon[i].1).powi(2))
                .sqrt();
        j = i;
    }
    perimeter
}

/**
Compute the equivalent perimeter of a polygon

```rust
use geometric_features::equivalent_perimeter;
let polygon: Vec<(f64,f64)> = vec![
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.0),
    (0.0, 1.0),
];
assert_eq!(equivalent_perimeter(&polygon), 3.5449077018110318);
```
 */
pub fn equivalent_perimeter(polygon: Polygon) -> f64 {
    (area(polygon) * consts::PI).sqrt() * 2.0
}

/**
Compute the compacity of a polygon

```rust
use geometric_features::compacity;
let polygon: Vec<(f64,f64)> = vec![
    (0.0, 0.0),
    (1.0, 0.0),
    (1.0, 1.0),
    (0.0, 1.0),
];
assert_eq!(compacity(&polygon), 0.10132118364233779);
```
 */
pub fn compacity(polygon: Polygon) -> f64 {
    area(polygon) / (perimeter(polygon).powi(2) / 4.0 * consts::PI)
}
