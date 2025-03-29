use std::io;
use std::any::type_name;

const max_number_of_equations:usize = 4;

struct LS {
    n: usize,
    A: [[f32;max_number_of_equations];max_number_of_equations],
    A_values : [f32;max_number_of_equations],       //It is assumed that a sparse matrix does not include non-zero values more that the number of rows or columns
    A_rows : [usize;max_number_of_equations],
    A_cols : [usize;max_number_of_equations],
    non_zero_elements_count: usize,
    b: [f32;max_number_of_equations],
    x: [f32;max_number_of_equations],
}

impl LS {
    fn modify_A_element(&mut self,row:usize,col:usize,new_value:f32) {
        self.A[row][col] = new_value;
    }
    fn modify_A_element_sparse(&mut self,row:usize,col:usize,new_value:f32) {
        self.A_values[self.non_zero_elements_count] = new_value;
        self.A_rows[self.non_zero_elements_count] = row;
        self.A_cols[self.non_zero_elements_count] = col;
        self.non_zero_elements_count += 1;
    }
    fn check_for_diagonal_dominant(&self) -> bool{
        let mut number_of_OK_rows = 0; // number of rows which the sum of off-diagonal elements is smaller than diagonal element.
        for row in 0..self.n {
            for col in 0..self.n {
                let mut sum0 = 0.0;
                for col in 0..self.n {
                    if row != col {
                        sum0 += self.A[row][col];
                    }
                }
                if sum0 < self.A[row][row] {
                    number_of_OK_rows += 1;
                }
            }
        }
        if number_of_OK_rows == self.n {
            return true;
        }
        else {
            return false;
        }
    }
    fn check_for_diagonal_dominant_sparse(&self) -> bool{
        let mut number_of_OK_rows = 0; // number of rows which the sum of off-diagonal elements is smaller than diagonal element.
        let mut sum_row = [0.0;max_number_of_equations];
        let mut diag_row = [0.0;max_number_of_equations];
        for i in 0..self.non_zero_elements_count {
            if self.A_rows[i] != self.A_cols[i] {
                sum_row[self.A_rows[i]] += self.A_values[i];
            } else {
                diag_row[self.A_rows[i]] = self.A_values[i];
            }
        }
        for i in 0..self.n {
            if sum_row[i] < diag_row[i] {
                number_of_OK_rows += 1;
            }
        }

        if number_of_OK_rows == self.n {
            return true;
        }
        else {
            return false;
        }
    }
    fn guass_sidel(&mut self) {
        const eps0:f32 = 0.0001;
        const max_itr:i16 = 1000;
        let mut err:f32 = 1.0;
        let mut itr: i16 = 0;
        
        while err > eps0 {
            for row in 0..self.n {
                let mut sum0 = 0.0;
                for col in 0..self.n {
                    if row != col {
                        sum0 += self.A[row][col] * self.x[col];
                    }
                }
                self.x[row] = (self.b[row]-sum0)/self.A[row][row];
            }
            err = 0.0;
            for row in 0..self.n {
                let mut sum0 = 0.0;
                for col in 0..self.n {
                    sum0 += self.A[row][col] * self.x[col];
                }
                err += (sum0 - self.b[row]).abs();
            }
            println!("itr = {} - err = {}",itr,err);
            itr += 1;        
            if itr > max_itr {
                break;
            }
        }
    }
    fn guass_sidel_sparse(&mut self) {
        const eps0:f32 = 0.0001;
        const max_itr:i16 = 1000;
        let mut err:f32 = 1.0;
        let mut itr: i16 = 0;

        while err > eps0 {
            let mut sum_row = [0.0;max_number_of_equations];
            let mut diag_row = [0.0;max_number_of_equations];
            for i in 0..self.non_zero_elements_count {
                if self.A_rows[i] != self.A_cols[i] {
                    sum_row[self.A_rows[i]] += self.A_values[i] * self.x[self.A_cols[i]];
                } else {
                    diag_row[self.A_rows[i]] = self.A_values[i] * self.x[self.A_cols[i]];
                }
            }
            for i in 0..self.n {
                self.x[i] = (self.b[i] - sum_row[i]) / diag_row[i];
            }
            let mut sum_row = [0.0;max_number_of_equations];
            let mut diag_row = [0.0;max_number_of_equations];
            for i in 0..self.non_zero_elements_count {
                sum_row[self.A_rows[i]] += self.A_values[i] * self.x[self.A_cols[i]];
            }                
            err = 0.0;
            for i in 0..self.non_zero_elements_count {
                err += (sum_row[i] - self.b[self.A_cols[i]]).abs();
            }
            println!("itr = {} - err = {}",itr,err);
            itr += 1;        
            if itr > max_itr {
                break;
            }
        }
    } 
    fn read_B(&mut self) {
        println!("Please enter the elements of constants vector (b). ");
        println!("The elements of a row could be seperated by comma.");
        let n = self.n;
        let mut input = String::new();
        io::stdin().read_line(&mut input).expect("Failed to get input");
        let input = input.trim();
        let mut col:usize = 0;
        for element in input.split(',') {
            let element_as_number: f32 = element.trim().parse().expect("Failed to convert to number");
            if col <= n {
                self.b[col] = element_as_number;
                col += 1;
            }
        }
        if col-1 >= n {
            println!("The maximum number of columns could not be more than {}. the elements after this item where trimed and neglected.",n);
        }
    }
    fn read_row_of_A(&mut self,row:usize) {
        let n = match row {
            0 => max_number_of_equations,
            other => self.n,
        };
        let mut input = String::new();
        io::stdin().read_line(&mut input).expect("Failed to get input");
        let input = input.trim();
        let mut col:usize = 0;
        for element in input.split(',') {
            let element_as_number: f32 = element.trim().parse().expect("Failed to convert to number");
            if col <= n {
                self.modify_A_element(row,col,element_as_number);
                // self.A[row][col] = element_as_number;
                col += 1;
            }
        }
        if col-1 >= n {
            println!("The maximum number of columns could not be more than {}. the elements after this item where trimed and neglected.",n);
        }
        if row == 0 {
            self.n = col;
        }
    }
    fn read_A_from_CLI(&mut self) {
        println!("Please enter the elements of coefficients matrix (A). For new row, press enter.");
        println!("The elements of a row could be seperated by comma.");
        self.read_row_of_A(0);
        for row in 1..self.n {
            self.read_row_of_A(row);
        }
    }
    
}




fn print_type_of<T>(_: &T) {
    println!("Type: {}", type_name::<T>());
}

fn main() {
    let mut ls0 = LS  {
        n: 2,
        A: [[2.0,3.0,0.0,0.0],[1.0,-1.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]],
        A_values: [0.0;max_number_of_equations],       //It is assumed that a sparse matrix does not include non-zero values more that the number of rows or columns
        A_rows: [0;max_number_of_equations],
        A_cols: [0;max_number_of_equations],
        non_zero_elements_count: 0,
        b: [5.0,3.0,0.0,0.0],
        x: [0.0,0.0,0.0,0.0],
    };
    ls0.read_A_from_CLI();
    if !ls0.check_for_diagonal_dominant() {
        println!("The Linear system is not diagonally dominant. so the gauess sidel may diverge.");
    }
    ls0.read_B();
    ls0.guass_sidel();
    println!("A = {:?}",ls0.A);
    println!("b = {:?}",ls0.b);
    println!("x = {:?}",ls0.x);
    println!("{}",ls0.n);

}



