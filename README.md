# Global Optimization
Repository for Global optimization learning work 

# Глобальная оптимизация
Репозиторий для хранения реализаций алгоритмов глобального поиска

## Задачи:

Реализовать:
* Последовательный алгоритм глобального поиска методом Стронгина для одномерных задач ( `OK` )
* Параллельный алгоритм глобального поиска методом Стронгина для одномерных задач ( `OK` )
* Параллельный алгоритм глобального поиска методом Стронгина для двумерных задач с использованием кривых Пеано (развертки) ( `OK` )
* Параллельный алгоритм глобального поиска методом Стронгина на двумерных и одномерных задач с применением техники GPGPU ( `OK` )
* * Тестирование АГП для двумерных задач на 100 функциях Гришагина ( `OK` )
* * Формульный транслятор ( `OK` )
* * Пользовательский интерфейс ( `OK` )
* * - Графопостроитель одномерных и двумерных функций ( `OK` )
* * - Формульный транслятор одномерных и двумерных функций ( `OK` )

## Описание репозитория:
* `AGS-sequental-one-size` -- последовательный АГП для одномерной задачи ( _+ параллельный по участкам интервала поиска_)
* `AGS-parallel-one-size` -- параллельный АГП для одномерной задачи ( _по характеристикам R для интервалов_ )
* `Optimizer` -- класс для хранения последовательного и параллельного АГП для одномерных, двумерных и N-мерных задач
* `Optimizer Visual` -- пользовательский оконный интерфейс прикладного характера в рамках интересов студентов
