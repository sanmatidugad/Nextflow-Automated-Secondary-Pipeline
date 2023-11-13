const csvParse = require('csv-parser');
const fs = require('fs');

// Main function to read and process the CSV file
(async (inputFilename, outputFilename) => {
    const parse = csvParse({});
    const genes = {};
    const data = [];
    let values = [];

    // Event listener for each row of data
    parse.on('data', (breaker) => {
        data.push(breaker);
    });

    // Event listener for the end of the file
    parse.on('end', async () => {
        // Capture the number of columns
        const nCols = Object.keys(data[0]).length;

        // Capture the original headers for the output file
        const headers = Object.keys(data[0]).join(',');

        // Process each column and row
        for (let col = 1; col < nCols; col++) {
            values.push([]);
            data.forEach((breaker, x) => {
                let columnKey = Object.keys(breaker)[col];
                let str = breaker[columnKey];

                const geneKeys = Object.keys(breaker).slice(0, 1);
                const gene = geneKeys.map(key => breaker[key]).join(',');

                const num = evenRound(parseFloat(str), 2);
                str = formatFloat(num);

                if (!genes[gene]) {
                    genes[gene] = [num];
                    values[values.length-1].push(str);
                } else if(genes[gene].length < values.length) {
                    genes[gene].push(num);
                    values[values.length-1].push(str);
                }
            });
        }

        // Sort the values
        values.forEach(L => {
            L.sort();
        });

        // Sort the genes
        let entries = Object.entries(genes);
        entries.sort((a, b) => (
            a[0].localeCompare(b[0])
        ));

        // Create a write stream for the output file
        const writeStream = fs.createWriteStream(outputFilename);

        // Write the original headers to the output file
        writeStream.write(headers + '\n');

        // Write the processed data to the output file
        for (let [value, L] of entries) {
            const TheCounts = [];
            for (let Y = 0; Y < L.length; Y++) {
                const norm = getNormalized(Y, parseInt(L[Y]), values);
                TheCounts[Y] = norm;  // Keep up to 3 decimal places
            }
            const CntStr = TheCounts.join(',');
            writeStream.write(value + ',' + CntStr + '\n');
        }
        writeStream.end();
    });

    // Read the input file and pipe it to the parser
    fs.createReadStream(inputFilename).pipe(parse);
})(process.argv[2], process.argv[3]);

// Function to normalize a value
const getNormalized = (col, n, values) => {
    const find = formatFloat(n);
    const i = values[col].indexOf(find);
    let sum = 0;

    if (i !== -1) {
        for (let x = 0; x < values.length; x++) {
            sum += parseInt(values[x][i]);
        }
    }

    // Keep up to 3 decimal places without rounding down
    return (sum / values.length).toFixed(3);
}

// Function to format a float as a string
const formatFloat = (num) => String(Math.round(num)).padStart(10, '0');

// Function to round a number to even or to a specified number of decimal places
const evenRound = (num, decimalPlaces = 0) => {
    const m = Math.pow(10, decimalPlaces);
    const n = +(decimalPlaces ? num * m : num).toFixed(12);
    const i = Math.floor(n), f = n - i;
    const e = 1e-8;
    const r = (f > 0.5 - e && f < 0.5 + e) ?
                ((i % 2 == 0) ? i : i + 1) : Math.round(n);
    return decimalPlaces ? r / m : r;
}
