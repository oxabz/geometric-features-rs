/*!
Contains the features related to the convex hull of a polygon
 */

use std::cmp::Ordering;

use crate::{area, perimeter, utils, Point, Polygon, PolygonOwned};

/**
Compute the convex hull of a polygon

The convex hull is the smallest convex polygon that contains all the points of the polygon
*/
fn convex_hull(polygon: Polygon) -> PolygonOwned {
    let mut convex_hull = Vec::with_capacity(polygon.len()); // The convex hull is at most the size of the polygon
    let mut polygon = polygon.to_vec(); // We need to sort the polygon, so we need to own it

    // Find the bottom point from which we start building the convex hull
    let min_i = utils::find_bottom(&polygon);
    let p = polygon.remove(min_i);
    convex_hull.push(p);

    // compute the relative point to the bottom point
    let relative_point = |pi: Point| -> (Point, Point) {
        let (x, y) = pi;
        (pi, (x - p.0, y - p.1))
    };

    // compute the angle between the bottom point and the other points
    let angle = |pi: (Point, Point)| -> (Point, f64, f64) {
        let length = (pi.1 .0.powi(2) + pi.1 .1.powi(2)).sqrt();
        let cos = pi.1 .0 / length;

        (pi.0, cos, length)
    };

    // Sort the points by the angle between the bottom point and the other points
    let sort = |(_, angle, length): &(Point, f64, f64),
                (_, oangle, olength): &(Point, f64, f64)|
     -> Ordering {
        match angle.partial_cmp(&oangle) {
            Some(Ordering::Equal) => olength.partial_cmp(&length),
            order => order,
        }
        .expect("Non orderable values!")
    };

    // Filter out the points that are on the same line as the bottom point
    let mut polygon: Vec<_> = polygon
        .into_iter()
        .filter(|&p| {
            let (x1, y1) = p;
            let (x2, y2) = convex_hull[0];
            (x1 - x2).abs() > 1e-6 || (y1 - y2).abs() > 1e-6
        })
        .map(relative_point)
        .map(angle)
        .collect();

    // Sort the polygon by the angle between the bottom point and the other points
    polygon.sort_by(sort);

    // Build the convex hull
    for (point, _, _) in polygon {
        // Backtrack until the angle is convex
        while convex_hull.len() > 1 {
            let (x1, y1) = convex_hull[convex_hull.len() - 2];
            let (x2, y2) = convex_hull[convex_hull.len() - 1];
            let (x3, y3) = point;

            let cross_product = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
            if cross_product > 0.0 {
                break;
            }
            convex_hull.pop();
        }
        convex_hull.push(point);
    }

    convex_hull.shrink_to_fit();
    convex_hull
}

/**
Compute the area of the convex hull of a polygon
 */
pub fn convex_area(polygon: Polygon) -> f64 {
    area(&convex_hull(polygon))
}

/**
Compute the perimeter of the convex hull of a polygon
 */
pub fn convex_perimeter(polygon: Polygon) -> f64 {
    perimeter(&convex_hull(polygon))
}

/**
Compute the deviation of the convex hull of a polygon
 */
pub fn convex_deviation(polygon: Polygon) -> f64 {
    let convex_hull = convex_hull(polygon);
    let convex_area = area(&convex_hull);
    let area = area(&convex_hull);
    (convex_area - area) / area
}

#[derive(Debug, PartialEq, Clone)]
pub struct ConvexHullFeatures {
    pub area: f64,
    pub perimeter: f64,
    pub deviation: f64,
}

/**
Compute the features of the convex hull of a polygon
 */
pub fn convex_hull_features(polygon: Polygon) -> ConvexHullFeatures {
    let convex_hull = convex_hull(polygon);
    println!("{:?}", convex_hull);
    let convex_area = area(&convex_hull);
    let area = area(&convex_hull);
    ConvexHullFeatures {
        area: convex_area,
        perimeter: perimeter(&convex_hull),
        deviation: (convex_area - area) / area,
    }
}
