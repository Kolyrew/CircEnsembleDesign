# Построчный разбор добавлений: кросс-пары IRES↔ORF (вероятностная структура)

Ниже — логическая последовательность: новые переменные и структуры → новые функции (построчно) → применение в коде.

---

## 1. Что такое threshold (порог)

**Threshold (порог)** — это число от 0 до 1 (например, 0.01), которое задаёт, **какие базовые пары считаются при подсчёте кросс-пар**.

- LinearPartition считает **матрицу вероятностей** P(i,j): для каждой пары позиций (i, j) — вероятность, что они спарены в равновесном ансамбле структур.
- Пары с очень малой вероятностью (шум) мы не хотим учитывать. Поэтому в вызов LinearPartition передаётся параметр **-c threshold**: выдавать только пары с **P(i,j) > threshold**.
- Итог: **кросс-парами** считаются только те пары (i ∈ IRES, j ∈ ORF или наоборот), у которых **P(i,j) > threshold**. При threshold = 0.01 учитываются пары с вероятностью выше 1%.

---

## 2. Новые переменные (что добавлено)

| Переменная | Файл | Тип | Значение по умолчанию | Смысл |
|------------|------|-----|------------------------|-------|
| `cross_pair_prob_threshold` | optimizer.h, ensemble_design.cpp, EnsembleDesign.py | double / float | 0.01 | Порог вероятности P(i,j) для детекции кросс-пар: учитываются только пары с P > этого значения. |
| `max_cross_pairs` | optimizer.h, ensemble_design.cpp, EnsembleDesign.py | int | 5 | Максимально допустимое число кросс-пар; последовательность с большим числом помечается как REJECTED. |

В Python они задаются аргументами: `--cross_pair_prob_threshold`, `--max_cross_pairs`.

---

## 3. optimizer.h — структура и объявления

### 3.1 Новые поля класса Optimizer (строки 137–138)

```cpp
double cross_pair_prob_threshold = 0.01;
int max_cross_pairs = 5;
```

- Добавлены два поля: порог вероятности и максимальное число кросс-пар.
- Используются при вызове `count_cross_pairs_probabilistic` и при выводе OK/REJECTED.

### 3.2 Сигнатура конструктора (строка 151)

```cpp
Optimizer(..., double cross_pair_prob_threshold_ = 0.01, int max_cross_pairs_ = 5);
```

- В конструктор добавлены два последних аргумента с значениями по умолчанию.

### 3.3 Структура CrossPairInfo (строки 168–172)

```cpp
struct CrossPairInfo {
    int count = 0;      // количество кросс-пар с P(i,j) > threshold
    double max_prob = 0.0;   // максимальная P среди этих кросс-пар
    double sum_prob = 0.0;   // сумма P по всем кросс-парам
};
```

- Результат подсчёта кросс-пар: сколько их, максимальная и суммарная вероятность.

### 3.4 Объявление функции (строка 174)

```cpp
CrossPairInfo count_cross_pairs_probabilistic(const std::string& seq, double prob_threshold = 0.01);
```

- Принимает последовательность и порог, возвращает `CrossPairInfo`.

### 3.5 Функция apply_cross_region_penalty (строки 233–239)

```cpp
inline ScoreType apply_cross_region_penalty(IndexType i_abs, IndexType j_abs) const {
    return 0.0;
}
```

- Раньше здесь добавлялся энергетический штраф за пары в граничной зоне. Сейчас **всегда 0** — penalty отключён, управление кросс-парами только через вероятностную детекцию и отбраковку.

---

## 4. optimizer.cpp — функция count_cross_pairs_probabilistic (строки 1302–1330)

**Назначение:** по полной последовательности (IRES+ORF) получить матрицу вероятностей пар от LinearPartition и подсчитать кросс-пары (IRES↔ORF) с P(i,j) > threshold.

Построчно:

| Строка | Код | Что делает |
|--------|-----|------------|
| 1302 | `Optimizer::CrossPairInfo Optimizer::count_cross_pairs_probabilistic(const std::string& seq, double prob_threshold)` | Объявление: последовательность и порог вероятности. |
| 1303 | `CrossPairInfo info;` | Создаётся структура результата (count=0, max_prob=0, sum_prob=0). |
| 1304 | `if (fixed_prefix_len == 0) return info;` | Если IRES нет (длина префикса 0), кросс-пар не бывает — возвращаем пустой результат. |
| 1306 | `std::string tmp_file = "/tmp/ensemble_bpp_" + std::to_string(getpid()) + ".txt";` | Имя временного файла для вывода LinearPartition (уникально по PID процесса). |
| 1307–1308 | `std::string cmd = "echo " + seq + " | ./tools/LinearPartition/linearpartition -V -c " + std::to_string(prob_threshold) + " -r " + tmp_file + " > /dev/null 2>&1";` | Строка команды: последовательность в stdin, `-V` (partition function), `-c prob_threshold` (только пары с P > threshold), `-r tmp_file` (запись пар в файл), вывод программы подавляется. |
| 1309 | `std::system(cmd.c_str());` | Запуск LinearPartition; в tmp_file оказываются строки вида "i j P(i,j)". |
| 1311 | `std::ifstream bpp_file(tmp_file);` | Открытие файла с матрицей вероятностей пар. |
| 1312 | `if (!bpp_file.is_open()) return info;` | Если файл не открылся (ошибка вызова), возвращаем пустой info. |
| 1314–1315 | `int i, j;` `double prob;` | Переменные для номера позиции i, j и вероятности. |
| 1316 | `while (bpp_file >> i >> j >> prob)` | Чтение каждой строки "i j prob" (позиции 1-based). |
| 1317–1318 | `bool i_in_ires = (i >= 1 && i <= fixed_prefix_len);` `bool j_in_ires = ...` | Проверка: лежит ли i (и j) в диапазоне IRES [1, fixed_prefix_len]. |
| 1319–1320 | `bool i_in_orf = (i > fixed_prefix_len);` `bool j_in_orf = ...` | Проверка: лежит ли позиция в ORF (после IRES). |
| 1321 | `if ((i_in_ires && j_in_orf) || (i_in_orf && j_in_ires))` | Пара (i,j) считается кросс-парой, если один конец в IRES, другой в ORF. |
| 1322 | `info.count++;` | Увеличиваем счётчик кросс-пар. |
| 1323 | `info.sum_prob += prob;` | Накапливаем сумму вероятностей. |
| 1324 | `if (prob > info.max_prob) info.max_prob = prob;` | Обновляем максимальную вероятность. |
| 1326 | `bpp_file.close();` | Закрываем файл. |
| 1327 | `std::remove(tmp_file.c_str());` | Удаляем временный файл. |
| 1329 | `return info;` | Возвращаем структуру с count, max_prob, sum_prob. |

**Итог:** функция по seq и threshold возвращает число кросс-пар, max P и sum P среди них.

---

## 5. optimizer.cpp — конструктор Optimizer (строки 1331–1335)

```cpp
Optimizer::Optimizer(..., double cross_pair_prob_threshold_, int max_cross_pairs_)
    : ... , cross_pair_prob_threshold(cross_pair_prob_threshold_), max_cross_pairs(max_cross_pairs_) {
    func1();
}
```

- Новые параметры сохраняются в поля `cross_pair_prob_threshold` и `max_cross_pairs`.

---

## 6. optimizer.cpp — цикл optimize(): мониторинг кросс-пар (строки 1232–1296)

### 6.1 Перед циклом (1232–1237)

| Строка | Код | Что делает |
|--------|-----|------------|
| 1232–1235 | `if (fixed_prefix_len > 0) { cout << "Cross-pair detection: prob_threshold = " << ... << ", max_cross_pairs = " << ... << endl; }` | Если есть IRES, выводим используемые порог и лимит кросс-пар. |
| 1237 | `int check_interval = 3;` | Проверять кросс-пары каждые 3 эпохи и на последней. |

### 6.2 Внутри цикла по эпохам (1266–1286)

| Строка | Код | Что делает |
|--------|-----|------------|
| 1266 | `if (fixed_prefix_len > 0 && (k % check_interval == 0 \|\| k == num_epochs))` | Условие: есть IRES и текущая эпоха — кратна 3 или последняя. |
| 1267–1268 | `auto current_seq = dfa.get_best_nuc_sequence(...);` `std::string full_seq = ires + current_seq.first;` | Берём текущую лучшую нуклеотидную последовательность по DFA и дописываем спереди IRES. |
| 1269 | `CrossPairInfo cp_info = count_cross_pairs_probabilistic(full_seq, cross_pair_prob_threshold);` | **Применение функции:** считаем кросс-пары по матрице вероятностей с заданным порогом. |
| 1271–1275 | `cout << "Iteration [" << ... << k << "]: obj: " << -end_state.inside << " \| cross-pairs(P>" << cross_pair_prob_threshold << "): " << cp_info.count << " max_P=" << cp_info.max_prob << " sum_P=" << cp_info.sum_prob` | Вывод номера итерации, целевой функции и статистики кросс-пар. |
| 1277–1281 | `if (cp_info.count > max_cross_pairs) { cout << " \| REJECTED"; } else { cout << " \| OK"; }` | Если кросс-пар больше допустимого — REJECTED, иначе OK (только вывод, параметры оптимизации не меняются). |
| 1283–1286 | `} else { cout << "Iteration [" << ... << k << "]: obj: " << -end_state.inside << endl; }` | В остальные эпохи выводим только номер итерации и obj. |

### 6.3 После цикла (1289–1296)

| Строка | Код | Что делает |
|--------|-----|------------|
| 1289–1292 | `if (fixed_prefix_len > 0) { ... final_seq_pair = dfa.get_best_nuc_sequence(...); final_full = ires + ...; final_cp = count_cross_pairs_probabilistic(final_full, ...);` | Снова **применение** count_cross_pairs_probabilistic к финальной лучшей последовательности. |
| 1293–1295 | `cout << "Final cross-pair analysis: " << final_cp.count << " pairs with P>" << ... << " (max_P=" << ... << ", sum_P=" << ... << ")" << endl;` | Вывод финальной статистики кросс-пар. |

**Итог:** внутри градиентного спуска каждые 3 эпохи и на последней вызывается `count_cross_pairs_probabilistic`, результат только выводится (OK/REJECTED), penalty не применяется.

---

## 7. ensemble_design.cpp — аргументы командной строки

| Строка | Код | Что делает |
|--------|-----|------------|
| 23–24 | `double cross_pair_prob_threshold = 0.01;` `int max_cross_pairs = 5;` | Локальные переменные со значениями по умолчанию. |
| 35–36 | `if (argc > 9) cross_pair_prob_threshold = atof(argv[9]);` `if (argc > 10) max_cross_pairs = atoi(argv[10]);` | Чтение 9-го и 10-го аргументов (после ires, ires_orf_lambda). |
| 77 | `Optimizer parser(..., cross_pair_prob_threshold, max_cross_pairs);` | Передача этих значений в конструктор Optimizer. |

---

## 8. EnsembleDesign.py — функция count_cross_pairs_probabilistic (строки 75–105)

**Назначение:** то же, что в C++ — по последовательности и порогу вернуть число кросс-пар, список пар с вероятностями и сумму вероятностей.

| Строка | Код | Что делает |
|--------|-----|------------|
| 75 | `def count_cross_pairs_probabilistic(seq, ires_len, prob_threshold=0.01):` | Имя функции, аргументы: seq, длина IRES, порог (по умолчанию 0.01). |
| 76 | `"""Count IRES-ORF cross-pairs using base pair probability matrix P(i,j) > threshold."""` | Документация. |
| 77–78 | `if ires_len == 0: return 0, [], 0.0` | Нет IRES — кросс-пар нет, возврат (count, список пар, sum_prob). |
| 79–81 | `import tempfile` `with tempfile.NamedTemporaryFile(..., delete=False) as tmp:` `tmp_path = tmp.name` | Создаётся временный файл для вывода LinearPartition. |
| 83 | `full_command = f"cd ./tools/LinearPartition && echo {seq} \| ./linearpartition -V -c {prob_threshold} -r {tmp_path}"` | Команда: LinearPartition с partition function (-V), порогом -c и записью пар в tmp_path. |
| 84 | `subprocess.run(full_command, shell=True, capture_output=True, text=True)` | Запуск команды. |
| 85–86 | `cross_pairs = []` `sum_prob = 0.0` | Список кросс-пар и сумма P. |
| 87–91 | `with open(tmp_path, 'r') as f:` `for line in f:` `parts = line.strip().split()` `if len(parts) != 3: continue` `i, j, prob = int(parts[0]), int(parts[1]), float(parts[2])` | Чтение файла построчно, разбор формата "i j prob". |
| 93–96 | `i_in_ires = (1 <= i <= ires_len)` ... `j_in_orf = (j > ires_len)` | Определение, в IRES или ORF лежат i и j. |
| 97 | `if (i_in_ires and j_in_orf) or (i_in_orf and j_in_ires):` | Условие кросс-пары. |
| 98–99 | `cross_pairs.append((i, j, prob))` `sum_prob += prob` | Добавление пары в список и накопление суммы. |
| 100 | `return len(cross_pairs), cross_pairs, sum_prob` | Возврат: количество, список, сумма. |
| 101–105 | `except Exception: return 0, [], 0.0` `finally: ... os.remove(tmp_path)` | При ошибке — пустой результат; временный файл всегда удаляется. |

---

## 9. EnsembleDesign.py — применение: команда запуска и отбор по кросс-парам

### 9.1 Формирование команды (строка 118)

```python
command = f"echo {protein} | {prog_path} ... {args.ires_orf_lambda} {args.cross_pair_prob_threshold} {args.max_cross_pairs}"
```

- В конец команды для C++ binary добавлены `cross_pair_prob_threshold` и `max_cross_pairs` (argv[9], argv[10]).

### 9.2 Аргументы парсера (строки 243–244)

```python
parser.add_argument('--max_cross_pairs', type=int, default=5, ...)
parser.add_argument('--cross_pair_prob_threshold', type=float, default=0.01, ...)
```

- Пользователь может задать порог и лимит кросс-пар из командной строки.

### 9.3 Отбор лучшей последовательности (строки 196–224)

| Строка | Код | Что делает |
|--------|-----|------------|
| 196 | `prob_thresh = args.cross_pair_prob_threshold` | Берём порог из аргументов. |
| 199–207 | Цикл по результатам каждого run: `cp_count, cp_list, sum_prob = count_cross_pairs_probabilistic(seq, ires_len, prob_thresh)` | **Применение** Python-версии функции к каждой финальной последовательности. |
| 205 | `status = "OK" if cp_count <= max_cp else "REJECTED"` | Статус по сравнению с max_cross_pairs. |
| 206 | `logs.append(... cross-pairs(P>{prob_thresh}): {cp_count} sum_P= ...)` | В лог пишется порог, число кросс-пар и sum_P. |
| 207 | `candidates.append((run_id, seq, efe, cp_count, cp_list))` | Кандидаты хранят в том числе cp_count. |
| 209 | `sorted_by_efe = sorted(candidates, key=lambda x: x[2])` | Сортировка по EFE (лучшая энергия сначала). |
| 211–215 | `for ... in sorted_by_efe:` `if cp_count <= max_cp:` `best_seq, best_efe, best_cp = seq, efe, cp_count` `break` | Выбор первого кандидата с числом кросс-пар не больше лимита. |
| 217–220 | `if best_seq is None ...: fallback = min(..., key=lambda x: (x[3], x[2]))` | Если все с превышением лимита — берём с минимальным числом кросс-пар, затем по EFE. |

**Итог:** после всех запусков для каждой последовательности вызывается `count_cross_pairs_probabilistic` с выбранным threshold; лучшая последовательность выбирается с учётом `max_cross_pairs` и EFE.

---

## 10. Краткая логическая цепочка

1. **Переменные:** вводятся `cross_pair_prob_threshold` (порог P) и `max_cross_pairs` (лимит числа кросс-пар); передаются из argv/argparse в Optimizer и в Python-логику отбора.
2. **CrossPairInfo** — структура результата: count, max_prob, sum_prob.
3. **count_cross_pairs_probabilistic (C++)** — по seq и threshold запускает LinearPartition, парсит файл с парами (i, j, P), отбирает кросс-пары (IRES↔ORF), считает count/max_prob/sum_prob.
4. **Применение в C++:** в optimize() каждые 3 эпохи и на последней вызывается count_cross_pairs_probabilistic для текущей лучшей последовательности; результат выводится (OK/REJECTED), penalty не используется.
5. **count_cross_pairs_probabilistic (Python)** — та же логика для финальных последовательностей каждого run.
6. **Применение в Python:** после всех run для каждого результата вызывается эта функция; по count и max_cross_pairs выбирается лучшая последовательность (с fallback, если все превышают лимит).
7. **apply_cross_region_penalty** — всегда возвращает 0; энергетический penalty отключён.

Таким образом, вся логика запрета кросс-парирования на текущий момент основана на **вероятностной структуре** (P(i,j) > threshold) и **отбраковке/выборе последовательностей**, без изменения энергии в DP.
