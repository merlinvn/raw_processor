use csv::{Reader, ReaderBuilder, Writer};
use serde::Deserialize;
use std::{error::Error, fs::File};

const PFPR_COL: usize = 12;
const NEW_INFECTIONS_COL: usize = 16;
const POSITIVE_CASES_COL: usize = 280;
const GENOTYPE0_COL: usize = 22;
const GENOTYPE127_COL: usize = 22 + 127;
const YEAR_COL: usize = 2;
const MONTH_COL: usize = 3;
const DAY_COL: usize = 4;

const C580Y_IDS: [usize; 64] = [
    4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23, 28, 29, 30, 31, 36, 37, 38, 39, 44, 45, 46, 47, 52,
    53, 54, 55, 60, 61, 62, 63, 68, 69, 70, 71, 76, 77, 78, 79, 84, 85, 86, 87, 92, 93, 94, 95,
    100, 101, 102, 103, 108, 109, 110, 111, 116, 117, 118, 119, 124, 125, 126, 127,
];

const PLAS_IDS: [usize; 64] = [
    2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 30, 31, 34, 35, 38, 39, 42, 43, 46, 47, 50,
    51, 54, 55, 58, 59, 62, 63, 66, 67, 70, 71, 74, 75, 78, 79, 82, 83, 86, 87, 90, 91, 94, 95, 98,
    99, 102, 103, 106, 107, 110, 111, 114, 115, 118, 119, 122, 123, 126, 127,
];

const KAF_OZ_IDS: [usize; 64] = [
    1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49,
    51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97,
    99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127,
];

const MDR2_IDS: [usize; 64] = [
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
    56, 57, 58, 59, 60, 61, 62, 63, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
    109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
];

#[derive(Debug, Deserialize, Eq, PartialEq)]
struct Row {
    values: Vec<String>,
}

struct SingleRunData {
    pfpr_at_year_10: Option<f64>,
    yearly_pfpr: Vec<(i32, f64)>,
}

fn extract_data(reader: &mut Reader<File>) -> Result<SingleRunData, Box<dyn Error>> {
    // for each row in the reader
    let iter = reader.deserialize();
    let mut run_result = SingleRunData {
        pfpr_at_year_10: None,
        yearly_pfpr: Vec::new(),
    };
    for result in iter {
        let record: Row = result.unwrap();
        let year = record.values[YEAR_COL].parse::<i32>().unwrap();
        let month = record.values[MONTH_COL].parse::<u32>().unwrap();
        let day = record.values[DAY_COL].parse::<u32>().unwrap();
        // check if year is 2032 month is 1 and day is 1
        if year == 2032 && month == 1 && day == 1 {
            run_result.pfpr_at_year_10 = Some(record.values[PFPR_COL].parse::<f64>().unwrap());
        }
        if month == 1 && day == 1 {
            run_result
                .yearly_pfpr
                .push((year - 2022, record.values[PFPR_COL].parse::<f64>().unwrap()));
        }
    }
    Ok(run_result)
}

fn main() {
    let mut pfpr_yearly = Writer::from_path("pfpr_yearly.csv").unwrap();

    // read folder path from command line
    let folder_path = std::env::args().nth(1).unwrap();
    let n_params = std::env::args().nth(2).unwrap().parse::<usize>().unwrap();
    let n_runs = std::env::args().nth(3).unwrap().parse::<usize>().unwrap();

    // let folder_path = "/home/neo/Projects/psu/MDA_AXX/A2_small/raw";

    for param_id in 0..n_params {
        for run_id in 0..n_runs {
            let job_id = param_id * 1000 + run_id;
            let file_path = format!("{}/monthly_data_{}.txt", folder_path, job_id);
            let file = File::open(file_path).expect("Unable to open file");
            let mut rdr = ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(file);

            let result = extract_data(&mut rdr);
            match result {
                Ok(result) => {
                    for (year, pfpr) in result.yearly_pfpr {
                        pfpr_yearly
                            .write_record(&[
                                param_id.to_string(),
                                run_id.to_string(),
                                year.to_string(),
                                pfpr.to_string(),
                            ])
                            .unwrap();
                    }
                }
                Err(e) => {
                    println!("param_id: {}, run_id:{}, Error: {}", param_id, run_id, e);
                    continue;
                }
            }
            // param id, run id (0-99), pfpr_at_year_10
            pfpr_yearly.flush().unwrap();
        }
    }
}
