/*!
Contains the feature function related to the related elipse of a polygon
 */

use nalgebra::{Matrix1xX, Matrix3, Matrix3xX};

use crate::{
    utils::{self},
    Point, Polygon,
};

const C_INV: [f64; 3 * 3] = [0.0, 0.0, 0.5, 0.0, -1.0, 0.0, 0.5, 0.0, 0.0];

/**
Fit an elipse to a polygon using a scanline algorythe and the "Numerically stable direct least squares fitting of ellipses" algorythm by Halir and Flusser

## Example

```rust
use geometric_features::related_ellipse::fit_ellipse;
let polygon: Vec<(f64,f64)> = vec![
    ( 8.04,  2.68),
    ( 5.4,   4.4 ),
    ( 6.1,   6.5 ),
    ( 9.24,  6.22),
    (12.12,  5.26),
    ( 9.86,  4.38),
    ( 10.0,  2.42),
    ( 8.86,  4.14)
];
let res = fit_ellipse(&polygon);
assert_eq!(res, (2.376780621008688, 2.2934319562351835, 3.999506964053637, (8.487782302188705, 5.068631718638504)));
```

 */
#[allow(non_snake_case)]
pub fn fit_ellipse(polygon: Polygon) -> (f64, f64, f64, Point) {
    let mut polygon = polygon.to_vec();
    polygon.push(polygon[0]);
    let interpolate = |arr:&[(f64, f64)]|{
        let (x0, y0) = arr[0];
        let (x1, y1) = arr[1];

        (0..5)
            .map(|x|x as f64 / 5.0)
            .map(move|t|(x0 * (1.0 - t) + x1 * t, y0 * (1.0 - t) + y1 * t))
    };
    let (x, y): (Vec<_>, Vec<_>) = polygon.windows(2).flat_map(interpolate).unzip();

    let x = Matrix1xX::from_row_slice(&x);
    let y = Matrix1xX::from_row_slice(&y);
    let n = x.len();

    let D1 = Matrix3xX::from_rows(&[
        x.component_mul(&x),
        x.component_mul(&y),
        y.component_mul(&y),
    ]);
    let D2 = Matrix3xX::from_rows(&[x, y, Matrix1xX::from_element(n, 1.0)]);

    let S1 = &D1 * D1.transpose();
    let S2 = D1 * D2.transpose();
    let S3 = &D2 * D2.transpose();

    let T = (-S3.try_inverse().unwrap()) * S2.transpose();
    let M = S1 + S2 * T;
    let Cinv = Matrix3::from_row_slice(&C_INV);
    let M = Cinv * M;

    let eigen = M.try_symmetric_eigen(1e-7, 200).unwrap();
    let (eigvec1, eigvec2, eigvec3) = (
        eigen.eigenvectors.column(0),
        eigen.eigenvectors.column(1),
        eigen.eigenvectors.column(2),
    );

    let con = 4.0
        * eigen
            .eigenvectors
            .row(0)
            .component_mul(&eigen.eigenvectors.row(2))
        - eigen
            .eigenvectors
            .row(1)
            .component_mul(&eigen.eigenvectors.row(1));
    let mut ak = vec![];
    if con[0] > 0.0 {
        ak.push(eigvec1);
    }
    if con[1] > 0.0 {
        ak.push(eigvec2);
    }
    if con[2] > 0.0 {
        ak.push(eigvec3);
    }
    let ak = Matrix3xX::from_columns(&ak);

    let akT = T * &ak;

    let a = ak[0];
    let b = ak[1] * 0.5;
    let c = ak[2];
    let d = akT[0] * 0.5;
    let f = akT[1] * 0.5;
    let g = akT[2];

    let den = b * b - a * c;

    assert!(den < 0.0, "Coefs don't match the elipse equation");

    let x0 = (c * d - b * f) / den;
    let y0 = (a * f - b * d) / den;

    let num = 2.0 * (a * f.powi(2) + c * d.powi(2) + g * b.powi(2) - 2.0 * b * d * f - a * c * g);
    let fac = ((a - c).powi(2) + 4.0 * b.powi(2)).sqrt();
    let a_axis = (num / (den * (fac - a - c))).sqrt();
    let b_axis = (num / (den * (-fac - a - c))).sqrt();

    let (width_gt_height, major_axis, minor_axis) = if a_axis > b_axis {
        (true, a_axis, b_axis)
    } else {
        (false, b_axis, a_axis)
    };

    let mut orientation = if b.abs() < 1e-6 {
        0.0
    } else if a <= c {
        (0.5 * (c - a + fac) / b).atan()
    } else {
        (0.5 * (c - a - fac) / b).atan() + std::f64::consts::PI / 2.0
    };

    if width_gt_height {
        orientation += std::f64::consts::PI / 2.0;
    }
    
    (major_axis, minor_axis, orientation, (x0, y0)) 
}

/**
Compute the eccentricity of the related elipse of a polygon

```rust
use geometric_features::{related_ellipse::eccentricity, Point, Polygon};
let polygon: Vec<(f64,f64)> = vec![
    ( 8.04,  2.68),
    ( 5.4,   4.4 ),
    ( 6.1,   6.5 ),
    ( 9.24,  6.22),
    (12.12,  5.26),
    ( 9.86,  4.38),
    ( 10.0,  2.42),
    ( 8.86,  4.14)
];
assert_eq!(eccentricity(&polygon), 0.26249954212114796);
```
 */
pub fn eccentricity(polygon: Polygon) -> f64 {
    let (major_axis, minor_axis, _, _) = fit_ellipse(polygon);
    (1.0 - minor_axis.powi(2) / major_axis.powi(2)).sqrt()
}

/**
Compute the minor and major axis of the related elipse of a polygon

```rust
use geometric_features::related_ellipse::principal_axes;

let polygon = [
    (-10.34, -11.32),
    (  3.61,  18.34),
    ( -3.07,  -5.34),
    ( 19.67,   6.5 )
];

let  (major_axis, minor_axis) = principal_axes(&polygon);
assert_eq!(major_axis, 23.247229367609776);
assert_eq!(minor_axis, 18.561527921248285);

```
*/
pub fn principal_axes(polygon: Polygon) -> (f64, f64) {
    let (major_axis, minor_axis, _, _) = fit_ellipse(polygon);
    (major_axis * 2.0, minor_axis * 2.0)
}

fn eliptic_deviation_(
    polygon: Polygon,
    precision: usize,
    hmajor_axis: f64,
    hminor_axis: f64,
    orientation: f64,
    center: (f64, f64),
) -> f64 {
    // Transform the polygon so that it is centered on the origin, axis aligned and scaled to fit in a unit circle
    let poly_centered = utils::transform::translate(&polygon, -center.0, -center.1);
    let poly_rotated = utils::transform::rotate(&poly_centered, -orientation);
    let poly_scaled = utils::transform::scale(&poly_rotated, 1.0 / hmajor_axis, 1.0 / hminor_axis);

    // Compute the vertical bounds of the polygon
    let max_y = poly_scaled
        .iter()
        .map(|&(_, y)| y)
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let min_y = poly_scaled
        .iter()
        .map(|&(_, y)| y)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    // Compute the deffect of the polygon
    let mut delta_acc = 0.0;
    for y in (0..precision).map(|x| x as f64 / precision as f64 * (max_y - min_y) + min_y) {
        let mut intersections = utils::find_horizontal_intersect(&poly_scaled, y);
        if y < 1.0 && y > -1.0 {
            let cos = y.asin().cos();
            intersections.push(cos);
            intersections.push(-cos);
        }

        intersections.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut it = intersections.chunks_exact(2);
        while let Some(&[x1, x2]) = it.next() {
            delta_acc = delta_acc + x2 - x1;
        }
    }

    let line_width = (max_y - min_y) / precision as f64;

    delta_acc * line_width * hmajor_axis * hminor_axis / crate::area(polygon)
}

/**
Compute the deviation beween the polygon and the related elipse using a scanline algorithm.
`precision` is the number of lines used to compute the deviation. The more lines, the more precise the result will be.

```rust
use geometric_features::related_ellipse::eliptic_deviation;
let polygon: Vec<(f64,f64)> = vec![
    ( 8.04,  2.68),
    ( 5.4,   4.4 ),
    ( 6.1,   6.5 ),
    ( 9.24,  6.22),
    (12.12,  5.26),
    ( 9.86,  4.38),
    ( 10.0,  2.42),
    ( 8.86,  4.14)
];
assert_eq!(eliptic_deviation(&polygon, 100), 0.5378961702988894);
```
 */
pub fn eliptic_deviation(polygon: Polygon, precision: usize) -> f64 {
    let (major_axis, minor_axis, orientation, center) = fit_ellipse(polygon);
    eliptic_deviation_(
        polygon,
        precision,
        major_axis,
        minor_axis,
        orientation,
        center,
    )
}

/**
Contains the features of the related elipse of a polygon
 */
#[derive(Debug, PartialEq, Clone)]
pub struct ElipseFeatures {
    pub major_axis: f64,
    pub minor_axis: f64,
    pub orientation: f64,
    pub center: (f64, f64),

    pub eccentricity: f64,
    pub eliptic_deviation: f64,
}

/**
Compute the features of the related elipse of a polygon. Usefull to avoid recomputing the related elipse multiple times.

```rust
use geometric_features::related_ellipse::{elipse_features, ElipseFeatures};

let polygon = [
    (-10.34, -11.32),
    (  3.61,  18.34),
    ( -3.07,  -5.34),
    ( 19.67,   6.5 )
];

let polygon2: Vec<(f64,f64)> = vec![
    ( 8.04,  2.68),
    ( 5.4,   4.4 ),
    ( 6.1,   6.5 ),
    ( 9.24,  6.22),
    (12.12,  5.26),
    ( 9.86,  4.38),
    ( 10.0,  2.42),
    ( 8.86,  4.14)
];

let features = elipse_features(&polygon, 200);
assert_eq!(features.major_axis, 23.247229367609776);
assert_eq!(features.minor_axis, 18.561527921248285);
assert_eq!(features.orientation, 1.1140533990057633);
assert_eq!(features.center, (3.816556476304908, 2.4698470796765917));
assert_eq!(features.eccentricity, 0.6020738096588001);
assert_eq!(features.eliptic_deviation, 3.1033783265455486);

let features = elipse_features(&polygon2, 200);
assert_eq!(
    ElipseFeatures { major_axis: 4.753561242017376, minor_axis: 4.586863912470367, orientation: 3.999506964053637, center: (8.487782302188705, 5.068631718638504), eccentricity: 0.26249954212114796, eliptic_deviation: 0.5379829218188382 },
    features
);

```
 */
pub fn elipse_features(polygon: Polygon, precision: usize) -> ElipseFeatures {
    let (major_axis, minor_axis, orientation, center) = fit_ellipse(polygon);
    let eccentricity = (1.0 - minor_axis.powi(2) / major_axis.powi(2)).sqrt();
    let eliptic_deviation = eliptic_deviation_(
        polygon,
        precision,
        major_axis,
        minor_axis,
        orientation,
        center,
    );
    let major_axis = major_axis * 2.0;
    let minor_axis = minor_axis * 2.0;
    ElipseFeatures {
        major_axis,
        minor_axis,
        orientation,
        center,
        eccentricity,
        eliptic_deviation,
    }
}

#[cfg(test)]
mod test {
    use crate::PolygonOwned;

    use super::*;
    use std::time;

    const TEST_POLYGON: [Point; 4] = [
        (-10.34, -11.32),
        (3.61, 18.34),
        (-3.07, -5.34),
        (19.67, 6.5),
    ];

    #[test]
    fn regression_test_fit_elipse() {
        let (major_axis, minor_axis, orientation, center) = fit_ellipse(&TEST_POLYGON);

        assert_eq!(major_axis, 11.623614683804888);
        assert_eq!(minor_axis, 9.280763960624142);
        assert_eq!(orientation.to_degrees(), 63.83055791523414);
        assert_eq!(center, (3.816556476304908, 2.4698470796765917));
    }

    const REPETITIONS: usize = 200;
    const POLYGON_SIZES: [usize; 9] = [4, 8, 16, 32, 64, 128, 256, 512, 1024];

    #[test]
    fn bench_equivalent_elipse() {
        println!("Benchmarking equivalent_elipse");
        println!("size, time");
        for &size in POLYGON_SIZES.iter() {
            let polygon = build_random_polygon(size);
            let start = time::Instant::now();
            for _ in 0..REPETITIONS {
                fit_ellipse(&polygon);
            }
            let elapsed = start.elapsed();
            println!("{}, {}", size, elapsed.as_millis());
        }
    }

    fn build_random_polygon(size: usize) -> PolygonOwned {
        (0..size)
            .map(|x| x as f64 / size as f64 * 2.0 * std::f64::consts::PI)
            .map(|angle| (angle, rand::random::<f64>() * 100.0))
            .map(|x| (x.0.cos() * x.1, x.0.sin() * x.1))
            .collect::<Vec<_>>()
    }

    #[test]
    fn find_minor_major_axis() {
        for _ in 0..100 {
            let polygon = build_random_polygon(100);
            let (major_axis, minor_axis, orientation, center_o_mass) = fit_ellipse(&polygon);
            assert!(major_axis > minor_axis);

            let mut best = (f64::INFINITY, usize::MAX);

            for s in 1..100 {
                let major_axis = major_axis * (s as f64 / 50.0);
                let minor_axis = minor_axis * (s as f64 / 50.0);
                let elliptic_deviation = eliptic_deviation_(
                    &polygon,
                    200,
                    major_axis,
                    minor_axis,
                    orientation,
                    center_o_mass,
                );

                if elliptic_deviation < best.0 {
                    best = (elliptic_deviation, s);
                }
            }

            println!("{}, {}, {}", major_axis, best.1 as f64 / 50.0, best.0);
        }
    }
}
