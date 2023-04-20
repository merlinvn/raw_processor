use csv::{Reader, ReaderBuilder, WriterBuilder};
use serde::Deserialize;
use std::{error::Error, fs::File};

const PFPR_COL: usize = 16;
const MONTHLY_NEW_INFECTIONS_COL: usize = 17;
const MONTHLY_POSITIVE_CASES_COL: usize = 158;
const CUMULATIVE_NTF_COL: usize = 23;
const CUMULATIVE_TF_COL: usize = 24;
const POPULATION_SIZE_COL: usize = 8;
const EIR_COL: usize = 10;
const GENOTYPE0_COL: usize = 29;
const GENOTYPE127_COL: usize = 29 + 127;
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
struct InputRow {
    values: Vec<String>,
}

#[derive(serde::Serialize)]
struct OutputRow {
    param_id: usize,
    run_id: usize,
    year: i32,
    pfpr: f64,
    c580y_freq: f64,
    plas_freq: f64,
    kaf_oz_freq: f64,
    mdr2_freq: f64,
    ntf: f64,
    tf: f64,
}

struct YearlyData {
    year: i32,
    pfpr: f64,
    c580y_freq: f64,
    plas_freq: f64,
    kaf_oz_freq: f64,
    mdr2_freq: f64,
    ntf: f64,
    tf: f64,
}

struct SingleRunData {
    yearly_data: Vec<YearlyData>,
    monthly_data: Vec<(i32, i32, f64, f64)>,
}

fn extract_data(reader: &mut Reader<File>) -> Result<SingleRunData, Box<dyn Error>> {
    // for each row in the reader
    let iter = reader.deserialize();
    let mut run_result = SingleRunData {
        yearly_data: Vec::new(),
        monthly_data: Vec::new(),
    };
    let mut is_complete = false;
    for result in iter {
        let record: InputRow = result?;
        let year = record.values[YEAR_COL].parse::<i32>()?;
        let month = record.values[MONTH_COL].parse::<u32>()?;
        let day = record.values[DAY_COL].parse::<u32>()?;

        if month == 1 && day == 1 {
            let pfpr = record.values[PFPR_COL].parse::<f64>()?;

            let mut sum = 0.0;
            for i in GENOTYPE0_COL..GENOTYPE127_COL + 1 {
                let genotype = record.values[i].parse::<f64>()?;
                sum += genotype;
            }

            let mut sum_c580y = 0.0;
            for i in C580Y_IDS.iter() {
                let genotype = record.values[*i + GENOTYPE0_COL].parse::<f64>()?;
                sum_c580y += genotype;
            }

            let mut sum_plas = 0.0;
            for i in PLAS_IDS.iter() {
                let genotype = record.values[*i + GENOTYPE0_COL].parse::<f64>()?;
                sum_plas += genotype;
            }

            let mut sum_kaf_oz = 0.0;
            for i in KAF_OZ_IDS.iter() {
                let genotype = record.values[*i + GENOTYPE0_COL].parse::<f64>()?;
                sum_kaf_oz += genotype;
            }

            let mut sum_mdr2 = 0.0;
            for i in MDR2_IDS.iter() {
                let genotype = record.values[*i + GENOTYPE0_COL].parse::<f64>()?;
                sum_mdr2 += genotype;
            }

            let pop_size = record.values[POPULATION_SIZE_COL].parse::<f64>()?;
            let mut ntf = record.values[CUMULATIVE_NTF_COL].parse::<f64>()?;
            let mut tf = record.values[CUMULATIVE_TF_COL].parse::<f64>()?;
            ntf *= 100.0 / pop_size;
            tf *= 100.0 / pop_size;

            // println!("{} {} {}", year - 2022, pfpr, sum_c580y / sum);

            run_result.yearly_data.push(YearlyData {
                year: year - 2022,
                pfpr,
                c580y_freq: sum_c580y / sum,
                plas_freq: sum_plas / sum,
                kaf_oz_freq: sum_kaf_oz / sum,
                mdr2_freq: sum_mdr2 / sum,
                ntf,
                tf,
            });
        }

        if year == 2027 && month == 1 && day == 1 {
            is_complete = true;
            break;
        }
    }

    if !is_complete {
        return Err("Run is not complete".into());
    }

    Ok(run_result)
}

fn main() {
    let mut yearly_data_file = WriterBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .from_path("yearly_data.csv")
        .unwrap();

    // read folder path from command line
    let folder_path = std::env::args().nth(1).unwrap();
    let n_params_from = std::env::args().nth(2).unwrap().parse::<usize>().unwrap();
    let n_params_to = std::env::args().nth(3).unwrap().parse::<usize>().unwrap();
    let n_runs = std::env::args().nth(4).unwrap().parse::<usize>().unwrap();

    // let folder_path = "/home/neo/Projects/psu/MDA_AXX/A2_small/raw";

    for param_id in n_params_from..n_params_to {
        println!("extracting data for param_id: {}", param_id);
        for run_id in 0..n_runs {
            let job_id = param_id * 1000 + run_id;
            let file_path = format!("{}/monthly_data_{}.txt", folder_path, job_id);
            let file = File::open(file_path);
            let file = match file {
                Ok(file) => file,
                Err(_) => {
                    println!(
                        "ERROR: Unable to open file for param_id: {}, run_id: {}, job_id: {}",
                        param_id, run_id, job_id
                    );
                    continue;
                }
            };

            let mut rdr = ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(file);

            let result = extract_data(&mut rdr);

            match result {
                Ok(result) => {
                    for yearly_data in result.yearly_data {
                        yearly_data_file
                            .serialize(OutputRow {
                                param_id,
                                run_id,
                                year: yearly_data.year,
                                pfpr: yearly_data.pfpr,
                                c580y_freq: yearly_data.c580y_freq,
                                plas_freq: yearly_data.plas_freq,
                                kaf_oz_freq: yearly_data.kaf_oz_freq,
                                mdr2_freq: yearly_data.mdr2_freq,
                                ntf: yearly_data.ntf,
                                tf: yearly_data.tf,
                            })
                            .unwrap();
                        // yearly_data_file
                        //     .write_record(&[
                        //         param_id.to_string(),
                        //         run_id.to_string(),
                        //         year.to_string(),
                        //         pfpr.to_string(),
                        //         c580y_freq.to_string(),
                        //     ])
                        //     .unwrap();
                    }
                }
                Err(e) => {
                    println!(
                        "ERROR: param_id: {}, run_id:{}, Error: {}",
                        param_id, run_id, e
                    );
                    continue;
                }
            }
            // param id, run id (0-99), pfpr_at_year_10
            yearly_data_file.flush().unwrap();
        }
    }
}
