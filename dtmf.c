#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define QMAX 700
#define MAXIN 500
#define NSYMB 16



int tabFreq[] = {697, 770, 852, 941, 1209, 1336, 1477, 1633}; // таблица частот
int tabFreq1[] = {1209, 1336, 1477, 1633, // массивы соответствий частот и кодируемых символов
                  1209, 1336, 1477, 1633,
                  1209, 1336, 1477, 1633,
                  1209, 1336, 1477, 1633};
int tabFreq2[] = {697, 697, 697, 697,
                  770, 770, 770, 770,
                  852, 852, 852, 852,
                  941, 941, 941, 941};

char symbols[] = {'1', '2', '3', 'A',
                  '4', '5', '6', 'B',
                  '7', '8', '9', 'C',
                  '*', '0', '#', 'D'};
double SymbolLenSec = 0.08; // длина буфера для сигнала одного символа в секундах
double PauseLenSec = 0.007; // длина буфера для сигнала паузы в секундах
double MinFreqAmplitude = 3.0; // минимальный порог для опознования амплитуды отдельной частоты
double MaxPauseDev = 0.32; // максимальное стандартное отклонение сигнала в буфере паузы
                           // для опознования паузы


int N; // количество отсчетов
int frequence; // частота дискретизации
double* inp; // массив отсчетов

char* output; // строка декодированных символов
int output_len; // их количество


typedef struct queues { // структура очереди для реализации буферов
    double qu[QMAX]; // содержимое очереди
    int last, first; // индекс последнего и первого элемента
} queue;

void init(queue *q) { // инициализация очереди
    q->first = 1;
    q->last = 0;
}

void insertx(queue *q, double x) { // вставка элемента в конец очереди
    if(q->last < QMAX-1) {
        q->last++;
        q->qu[q->last]=x;
    }
    else
        printf("Очередь полна!\n");
}

int isempty(queue *q) { // проверка очереди на наличие элементов
    if(q->last < q->first)    return 1;
    else  return 0;
}

double removex(queue *q) { // удаление элемента с начала очереди
    double x;
    if(isempty(q)==1) {
        printf("Очередь пуста!\n");
        return 0;
    }
    x = q->qu[q->first];
    for(int h = q->first; h < q->last; h++) {
        q->qu[h] = q->qu[h+1];
    }
    q->last--;
    return x;
}

int len(queue *q) { // длина очереди
    return q->last - q->first + 1;
}

double goertzel(const queue *sample, int sample_len, int sample_freq, int target_freq) { // реализация алгоритма Герцеля
    int k;
    // номер частоты, для которой вычисляется амплитуда в данной выборке
    k = (int)((double)target_freq / sample_freq * sample_len) + 1;
    double s[2] = {0}; // массив вычислеямых промежуточных значений для алгоритма
    double t;
    double super_k = 2 * M_PI * k / sample_len; // просто постоянный коэффициент используемый в вычислениях
    for (int i = 1; i < sample_len + 1; i++) {
        t = 2 * cos(super_k) * s[1] - s[0] + sample->qu[i];
        s[0] = s[1];
        s[1] = t;
    }
    double t1 = s[1];
    double t2 = s[0];
    double amplitude = sqrt(t1 * t1 - 2 * cos(super_k) * t1 * t2 + t2 * t2);
    return amplitude;
}

void input_processing(char* text){ // обработка входных данных
    FILE *f;
    f = fopen(text, "r");
    if (f == NULL){
        printf("File reading error");
        return;
    }
    // заполняет переменные: массив отсчетов, их количество, и частоту дискретизации
    fscanf(f, "%d %d", &N, &frequence);
    inp = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++){
        fscanf(f, "%lf", &inp[i]);
    };
    fclose(f);
}

int pause_check(queue *pause, int pause_len){ // проверка буфера паузы для ее обнаружения
    double deviation = 0, s1 = 0, s2 = 0, mean, t;
    // считается стандартное отклонение отсчетов
    // и сравнивается с пороговым значением
    for (int i = 1; i < pause_len + 1; i++) {
        t = pause->qu[i];
        s1 += t;
        s2 += t*t;
    }
    mean = s1 / pause_len;
    deviation = sqrt((s2 + mean * mean * pause_len - 2 * s1 * mean) / pause_len);
    return deviation <= MaxPauseDev;
}


// функция для определения индекса максимального элемента в массиве double
int max_index(const double* arr, const int len){
    double max = arr[0];
    int max_ind = 0;
    for (int i = 0; i < len; i++){
        if (arr[i] > max) {
            max = arr[i];
            max_ind = i;
        }
    }
    return max_ind;
}



void make_symbol(double* comp_arr){ // функция для добавления символа
    int maxI = max_index(comp_arr, NSYMB); //  индекс символа с наибольшим счетчиком
    if (comp_arr[maxI] > 0){
        output = (char*)realloc(output, (output_len + 1) * sizeof(char));
        output[output_len] = symbols[maxI];
        output_len += 1;
        // очищаем словарь для поиска символов, т.к. сигнал, кодирующий символ завершен,
        // либо только начинается
        for (int i = 0; i < NSYMB; i++){
            comp_arr[i] = 0;
        }
    }
}

void find_amp(double* amps, queue* symb, int num, int freq){
    // нахождение амплитуд восьми нужных частот в данной выборке
    // (буфере символа)
    for (int i = 0; i < NSYMB/2; i++){
        amps[i] = goertzel(symb, num, freq, tabFreq[i]);
    }
}

void find_symb(double* amps, double* found_symb){
    // фунция находит две частоты с максимальными амплитудами и увеличивает
    // счетчик символа, соответствующего этим амплитудам
    int max1_ind = max_index(amps, NSYMB / 2);
    double max1 = amps[max1_ind];
    amps[max1_ind] = -1;
    int max2_ind = max_index(amps, NSYMB / 2);
    double max2 = amps[max2_ind];
    amps[max2_ind] = -1;
    int max3_ind = max_index(amps, NSYMB / 2);
    // проверка на шум, если есть третья "большая" частота, то это сокрей всего он
    if (amps[max3_ind] > max2 / 2){
        return;
    }
    double freq1, freq2;
    // также проверяем амплитуды частот
    if (max1 >= MinFreqAmplitude && max2 >= MinFreqAmplitude){
        if (max1_ind > max2_ind){
            freq1 = tabFreq[max1_ind];
            freq2 = tabFreq[max2_ind];
        }
        else{
            freq1 = tabFreq[max2_ind];
            freq2 = tabFreq[max1_ind];
        }
        for (int i = 0; i < NSYMB; i++){
            if (freq2 == tabFreq2[i] && freq1 == tabFreq1[i]){
                found_symb[i] += 1;
            }
        }
    }
}

void data_processing(double* data, int num, int freq) { // основная функция обработки массива данных
    // длины буферов в отсчетах
    int pause_len = (int)(PauseLenSec * freq);
    int symbol_len = (int)(SymbolLenSec * freq);
    // сами буфера
    queue *pause;
    queue *symbol;
    // массив, содержащий при обнаружении соответствующего символа счетчик его
    // обнаружения (символ определяется по индексу в массиве)
    double found_symbols[NSYMB] = {0};
    //  инициализация буферов, строки на вывод, массива амплитуд (меняется каждый раз,
    //  когда буфер символа заполняется)
    output = NULL;
    output_len = 0;
    output = (char*)realloc(output, (output_len + 1) * sizeof(char));
    double amplitudes[NSYMB/2] = {0};
    pause = (queue*)malloc(sizeof(queue));
    symbol = (queue*)malloc(sizeof(queue));
    init(pause);
    init(symbol);

    for (int i = 0; i < num; i++) {
        // сначала заполняется буфер паузы
        if (len(pause) < pause_len) {
            insertx(pause, data[i]);
        }
        else if (pause_check(pause, pause_len)) { //  если находится в паузе (тишине) или на переходе с тишины на сигнал символа
            init(pause); //  очищаются буфера, т. к. пауза найдена
            init(symbol); //  предыдущий блок отсчетов для детектирования символа уже был обработан либо не актуален
            make_symbol(found_symbols); //  добавляем задетектированный символ в список для вывода
            insertx(pause, data[i]); //  не забыть добавить новый отсчет в буфер для паузы
        }
        else{ // не пауза
            //  перекидываем четверть старых отсчетов из буфера паузы в буфер символа
            while (len(pause) > (int)(pause_len * 3 / 4) && len(symbol) < symbol_len){
                insertx(symbol, pause->qu[pause->first]);
                removex(pause);
            }
            if (len(symbol) >= symbol_len){ //  если длины буфера символа достаточно для его обнаружения,
                find_amp(amplitudes, symbol, len(symbol), freq); //  то вычисляем амплитуды для восьми необходимых частот
                find_symb(amplitudes, found_symbols); //  находится наиболее возможный символ
                while (len(symbol) > (int)(symbol_len * 3 / 4)){ //  удаляем четверть самых старых отсчетов
                    removex(symbol);
                }
            }
            insertx(pause, data[i]);
        }
    }
}


int main() {
    char input[MAXIN];
    printf("Enter the path to the decryption file:\n");
    gets(input);
    input_processing(input);
    data_processing(inp, N, frequence);
    printf("Decoding result:\n");
    for (int i = 0; i < output_len; i++){
        printf("%c", output[i]);
    }
    return 0;
}